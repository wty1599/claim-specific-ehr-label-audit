# D3 Model And Figure Source Scope

Date: 2026-10-08. Source lineage, not an independent clinical rerun.

| Model/source | Current source | Temperature-correction coverage |
| --- | --- | --- |
| Base | Retained archived comparator; no temperature input | Not newly fitted by temperature-correction or Raw-completion scripts. |
| Base + K2 | Archived fitted coefficients applied to corrected external L1 labels | Correction 03 updates predictions; 04-05 bootstrap/recalibrate. |
| Raw-EN and Raw-RF | Verified reconstructed source fits, Base covariates plus 33 features, unchanged archived L2 recipe | Raw completion 01 applies corrected temperature for both external outcomes and computes patient-bootstrap AUC and current paired contrasts. |
| Figure 4C | Corrected Base + K2 and corrected Raw results, with retained Base and source-internal CV estimates | Completed Raw panel-C delivery replaces earlier historical Raw rows; not four newly fitted external models. |
| Figure 4D | Corrected Base + K2 and retained Base, held-out mortality recalibration | Raw completion does not refit or alter this source. |

## Historical Persistence And Verified Reconstruction

The historical generator's `fit_store` was in memory. The historical
head-to-head script saved summary/preprocessing/partition provenance rather
than fitted outcome models. The completion therefore reconstructs models,
not retrieves an original archived Raw fit.

The authorized source run preserved the four-model/two-outcome order,
original seed, ten-fold loops and intervening external bootstrap calls.
Its reported patient-aligned historical prediction check passed within 1e-12
(maximum absolute error approximately 6.7e-16). Current Raw AUC and paired
contrasts each use 1,000 patient (`uniquepid`) bootstrap replicates per
statistic, with zero reported failures. These are execution-source records;
public packaging did not repeat them or inspect patient rows.

Raw preprocessing remains the archived L2 `object$recipe`. It must not be
replaced by the recovered L1 label-assignment coordinates. Both Raw models
use Base covariates plus the 33 features, including temperature. The current
paired source is `raw_vs_corrected_base_k2_paired_auc.csv` from Raw completion,
not the earlier correction-04 comparison with historical Raw probabilities.

## Calibration And Historical Comparison

The displayed intercept remains offset-only calibration-in-the-large (CITL);
the slope comes from a joint logistic fit. The separately named joint-fit
intercept does not replace Figure 4's convention. The dated addendum clarifies
the original specification without adding intervals or changing standards.

Completion 03 resolves tiny historical RF AUC rank-tie changes caused by
archived prediction-CSV rounding. It uses fully reproduced historical point
estimates and original intervals without refitting or repeating bootstrap.
The validated historical/completed comparison files are authoritative;
first-pass CSV-readback summaries remain diagnostics.

L1 assignment and Raw L2 preprocessing are distinct, explicitly retained
rules. Full-eICU and held-out-hospital calibration are different evaluation
populations. Current Figure 4A/B/D and plotting calculations remain unchanged.
