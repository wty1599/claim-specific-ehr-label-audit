# Claim-specific EHR label audit: research code

Code for *What a reused EHR-derived patient label can support: a claim-specific specification illustrated with a sepsis partition*.

Version `v1.0.0-jamia-submission` (2026-10-08).
Repository: <https://github.com/wty1599/claim-specific-ehr-label-audit>.

## Scientific Scope

The claim-specific specification is illustrated with one fixed sepsis K=2 partition. It is not a validation across multiple partitions or diseases. Stable assignment alone does not establish discrete biological classes, added predictive information, external semantic equivalence, or treatment-predictive value.

The current D3 source includes the post hoc 2026-10-08 eICU temperature correction. Fixed MIMIC-IV labels, recovered 33-feature coordinates and centroids remain unchanged. Valid nurse-charted Celsius values take priority, followed by converted Fahrenheit and valid periodic measurements; other external features remain unchanged.

Script 03 in `source_correction/` updates Base + K2 external predictions; Base has no temperature input and is retained. The subsequent authorized `raw_temperature_completion/` reconstructs the original source fits, verifies archived predictions for four models and two outcomes, and evaluates Raw-EN/Raw-RF with corrected temperature. Raw models use the unchanged archived L2 recipe, not the recovered L1 assignment scale. Current Figure 4C and current paired comparisons use the completed corrected Raw results, not historical Raw reference rows. No claim of four newly fitted corrected external models is made.

Model recovery/evaluation uses R 4.6.0; the unchanged plotting calculations use R 4.5.3. The original Raw fit objects were not saved, so the completed models are verified reconstructions rather than retrieved original objects. See `docs/D3_CALIBRATION_SOURCE_NOTE.md`.

## Code-Only Boundary

No clinical source records, patient-level derivatives, fitted patient-derived objects, aggregate result tables, hospital rows, figure datasets, rendered figures or other binary artifacts are distributed. MIT covers code/documentation only, not restricted data or derivatives. MIMIC-IV and eICU-CRD access and use conditions continue to apply.

Do not place restricted inputs or generated artifacts inside this directory. `config/paths.example.R` lists configuration variables. Correction scripts require existing external `EHR_AUDIT_RESTRICTED_DATA_ROOT` and `EHR_AUDIT_OUTPUT_ROOT` directories. `EHR_AUDIT_REPO_ROOT` selects this code directory.

## Code Map

| Path | Purpose |
| --- | --- |
| `sql/`, `R/partition/` | Extraction and fixed K=2 construction; historical eICU extraction precedes the documented temperature correction. |
| `R/step0/` | Clinical and semantic composition. |
| `R/domain1/`, `R/domain2/` | Discreteness controls and incremental-information comparisons. |
| `R/domain3/` | Historical comparator, parameter recovery, L1 deployment and current source correction. |
| `R/domain3/source_correction/` | Corrected assignment, historical Base + K2-only bootstrap/comparison stage, recalibration, funnel and missing-pattern transfer. |
| `R/domain3/raw_temperature_completion/` | Verified source-model reconstruction, corrected Raw evaluation, current paired comparisons and final Figure 4C source staging. |
| `sql/source_audit/` | Read-only reconciliation queries and restricted temperature export. |
| `R/domain4/`, `R/simulation/` | Support diagnostics and retained known-truth controls, not new treatment-effect analyses. |
| `figures/` | Main/supplementary renderers; one current Figure 4 renderer. |
| `tests/` | Static boundary checks, R parsing and small artificial-fixture tests. |

Main figure sources use `EHR_AUDIT_MAIN_FIGURE_SOURCE_ROOT`, defaulting to the external output root's `main_figures/source_data/`. Final Figure 4 uses `EHR_AUDIT_FIGURE4_SOURCE_ROOT`, defaulting to `domain3_raw_temperature_completion/figure4/source_data/`. An optional Figure 4 directory argument takes precedence in both `figures/main/scripts/run_all.R` and `figures/D3/scripts/Figure4.R`. Correction-stage preparation uses the separate `EHR_AUDIT_CORRECTION_FIGURE_SOURCE_ROOT`. Generated figures, QA and actual consumed-input fingerprints remain outside the repository.

See `docs/CLAIM_SCRIPT_MAP.md`, `docs/FIGURE_CROSSWALK.md` and `docs/SOURCE_PATH_CROSSWALK.md`. Excluded intermediate artifacts prevent one-command replay.

## Reproducibility Limits

No database query, analysis campaign, clinical model fit, bootstrap campaign or figure render was executed while preparing this directory. Earlier authorized calculations differ from this release's static checks. Original D2 ten-fold assignments were not saved, so exact fold-level replay is unavailable. Historical D1/simulation scripts require excluded intermediate inputs. Public path/language adaptations do not constitute fresh numerical validation.

Run `tests/check_code_release.ps1`, `Rscript tests/parse_release_scripts.R`, and `Rscript tests/test_public_dependencies.R` for limited local checks. Audit records remain outside the repository.

## License And Citation

MIT (`LICENSE`). `CITATION.cff` and `.zenodo.json` record confirmed author order and author-verified title-page affiliations, without ORCIDs or private email addresses. Cite the exact commit/archive version used.
