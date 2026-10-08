param(
    [ValidateSet('all', 'preflight', 'recover', 'evaluate', 'finalize', 'verify', 'stage-figure')]
    [string]$Stage = 'all',
    [string]$DataRoot = $env:EHR_AUDIT_RESTRICTED_DATA_ROOT,
    [string]$Rscript = $env:EHR_AUDIT_MODEL_RSCRIPT
)
$ErrorActionPreference = 'Stop'
$code = $PSScriptRoot
$repo = (Resolve-Path -LiteralPath (Join-Path $code '../../..')).Path
if (!$Rscript) { $Rscript = (Get-Command Rscript.exe -ErrorAction Stop).Source }
if (!(Test-Path -LiteralPath $Rscript)) { throw 'Required R runtime is unavailable.' }
foreach ($path in @($DataRoot, $env:EHR_AUDIT_OUTPUT_ROOT)) {
    if (!$path -or !(Test-Path -LiteralPath $path -PathType Container)) { throw 'Set existing external data and output roots.' }
    $resolved = (Resolve-Path -LiteralPath $path).Path.TrimEnd('\')
    if ($resolved -eq $repo -or $resolved.StartsWith($repo + '\', [StringComparison]::OrdinalIgnoreCase)) {
        throw 'Restricted inputs and generated outputs must remain outside the code repository.'
    }
}
$root = Join-Path $env:EHR_AUDIT_OUTPUT_ROOT 'domain3_raw_temperature_completion'
$logs = Join-Path $root 'logs'
New-Item -ItemType Directory -Path $logs -Force | Out-Null
$stages = if ($Stage -eq 'all') { @('preflight', 'recover', 'evaluate', 'finalize', 'verify', 'stage-figure') } else { @($Stage) }
$previous = $env:EHR_AUDIT_RESTRICTED_DATA_ROOT
$previousRepo = $env:EHR_AUDIT_REPO_ROOT
try {
    $env:EHR_AUDIT_RESTRICTED_DATA_ROOT = $DataRoot
    $env:EHR_AUDIT_REPO_ROOT = $repo
    foreach ($part in $stages) {
        $script = switch ($part) {
            'finalize' { '03_finalize_historical_comparison.R' }
            'verify' { '02_verify_completion.R' }
            'stage-figure' { '04_stage_current_figure4.R' }
            default { '01_restore_and_evaluate_raw.R' }
        }
        $entry = Join-Path $code ('scripts/' + $script)
        $arguments = if ($part -in @('preflight', 'recover', 'evaluate')) { @('--stage=' + $part) } else { @() }
        $log = Join-Path $logs ($part + '.log')
        if (Test-Path -LiteralPath $log) { throw "Refusing to overwrite run log: $log" }
        $ErrorActionPreference = 'Continue'
        & $Rscript --vanilla $entry @arguments 2>&1 | Tee-Object -FilePath $log
        $exitCode = $LASTEXITCODE
        $ErrorActionPreference = 'Stop'
        if ($exitCode -ne 0) { throw "R stage failed: $part (exit $exitCode)" }
    }
} finally {
    $env:EHR_AUDIT_RESTRICTED_DATA_ROOT = $previous
    $env:EHR_AUDIT_REPO_ROOT = $previousRepo
}
