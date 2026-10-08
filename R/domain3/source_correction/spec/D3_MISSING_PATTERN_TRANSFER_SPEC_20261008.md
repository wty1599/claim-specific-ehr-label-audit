# Post hoc eICU missing-pattern transfer check

Written 2026-10-08 after the original D3, recalibration, and temperature-source correction results were known. This is a sensitivity analysis of technical assignment, not prospective registration or a test of external semantic validity. The specification is frozen before computing its outcomes; SHA-256 is recorded separately.

## Inputs and frozen objects

- All 4,689 MIMIC-IV records with finite values for all 33 clustering inputs in the existing `final_full.csv` extract. Use their archived fixed K2 labels, not any new clustering.
- The 17,465 eICU stays' 33-input missingness masks after the specified nurse-charted temperature correction. A mask is nonfinite or missing in the corrected feature extract. Do not use eICU outcomes or labels to select masks.
- The already recovered original centroids, winsorization bounds, standardization, and source-median imputation values. Verify the exact parameter-file hash before analysis.

## Simulation

Use 500 independent repetitions, seed 20261008. In each repetition, sample 4,689 eICU stay-level missingness masks with replacement from all 17,465 corrected eICU masks and apply one sampled mask to each complete MIMIC-IV record. Preserve the MIMIC-IV observed values for unmasked inputs. Impute newly masked inputs with the fixed source medians, apply the recovered original coordinate system and nearest-centroid assignment, and compare with each record's archived fixed label. Do not rerun MICE, estimate centroids, or use outcomes. A complete-record baseline must reproduce all 4,689 archived labels before repetitions.

For each repetition report: total and directional C1-to-C2/C2-to-C1 flips; original and induced C1 prevalence; induced prevalence minus the same 4,689-record baseline; and prevalence/flip direction by assigned mask burden (0, 1-2, 3-5, at least 6 missing inputs). Retain all repetitions, including zero-flip ones. Aggregate across repetitions as mean, Monte Carlo standard error, and empirical 2.5th/97.5th percentiles. Compare the induced missingness-stratum prevalence gradient descriptively with the observed corrected eICU gradient. Do not fit outcome models or call the induced labels external accuracy.

## Interpretation boundary

Independent reassignment of eICU masks to complete MIMIC-IV records estimates mechanical label movement under the external pattern distribution. It does not preserve the joint distribution of eICU physiology and missingness, assess external case mix, or identify the fraction of actual eICU prevalence differences caused by missingness. The complete-record subset may differ substantially from the full MIMIC-IV cohort; use its own initial C1 prevalence as the reference. Patient-level masked values, labels, and sampled mask indices stay in the restricted workspace. Only aggregate summaries may be released after disclosure review.
