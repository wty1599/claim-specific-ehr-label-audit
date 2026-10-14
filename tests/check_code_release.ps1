param(
  [string]$RepositoryRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
)

$ErrorActionPreference = 'Stop'
$root = (Resolve-Path -LiteralPath $RepositoryRoot).Path
$issues = [System.Collections.Generic.List[string]]::new()
$allowedCsv = @(
  'R_PACKAGE_REQUIREMENTS.csv',
  'R/simulation/spec/D3_SIM_PUBLIC_SOURCE_HASHES_v1.csv',
  'R/simulation/spec/D4_SIM_PUBLIC_SOURCE_HASHES_v2.csv'
)
$allowedExtensions = @('.r', '.sql', '.md', '.csv', '.json', '.ps1', '.cff')
$allowedNames = @('LICENSE', '.gitignore', '.gitattributes')
$personalPath = '(?i)([a-z]:[/\\]users[/\\]|/users/[^/\s]+/|/home/[^/\s]+/)'
$credential = '(?i)(gh[pousr]_[a-z0-9]{30,}|github_pat_[a-z0-9_]{30,}|-----BEGIN (RSA |EC |OPENSSH )?PRIVATE KEY-----)'

$files = Get-ChildItem -LiteralPath $root -Recurse -File | Where-Object {
  -not $_.FullName.Substring($root.Length + 1).Replace('\', '/').StartsWith('.git/')
}
foreach ($file in $files) {
  $relative = $file.FullName.Substring($root.Length + 1).Replace('\', '/')
  if ($file.Attributes -band [IO.FileAttributes]::ReparsePoint) {
    $issues.Add("symlink: $relative")
  }
  if ($relative -match '^(data/|figure_projects/main/source_data/|figure_projects/supplement/data/)' -or
      $relative -match '(^|/)(data_restricted|patient_level|private_local|private_not_for_release|source_data|aggregate_outputs|generated|outputs|output|qa|logs|checkpoints|code_release_audit)/' -or
      $relative -match '(^|/)output_[^/]+/') {
    $issues.Add("restricted directory: $relative")
  }
  if ($file.Extension.ToLowerInvariant() -notin $allowedExtensions -and
      $file.Name -notin $allowedNames) {
    $issues.Add("non-code or unexpected extension: $relative")
  }
  if ($file.Extension.ToLowerInvariant() -eq '.csv' -and $relative -notin $allowedCsv) {
    $issues.Add("unexpected CSV: $relative")
  }
  if ($relative -eq 'tests/check_code_release.ps1') { continue }
  $content = Get-Content -LiteralPath $file.FullName -Raw -Encoding UTF8
  if ($content -match $personalPath) {
    $issues.Add("workstation path: $relative")
  }
  if ($content -match $credential) {
    $issues.Add("possible credential: $relative")
  }
  if ($content -match '[\p{IsCJKUnifiedIdeographs}]') {
    $issues.Add("Chinese text: $relative")
  }
  if ($content -match '(?i)(?:^|["''\s(])[a-z]:[/\\]|/home/|/Users/') {
    $issues.Add("absolute workstation path: $relative")
  }
  $bytes = [IO.File]::ReadAllBytes($file.FullName)
  if ($bytes -contains 0) { $issues.Add("binary content: $relative") }
  if ($bytes -contains 13) { $issues.Add("non-LF newline: $relative") }
}

if ($issues.Count) {
  $issues | Sort-Object -Unique | ForEach-Object { Write-Error $_ }
  exit 1
}
Write-Host "Code-release scan passed for $($files.Count) files. This is not a privacy certification."
