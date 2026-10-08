# D3 analysis generations

The archived `generator/` analysis is the historical reconstructed L2 rule.
The newer `parameter_recovery/` scripts recover the original fixed-label
centroids and coordinate parameters without reclustering; `deployed_rule/`
uses those parameters with fixed source-median imputation (L1); and
`recalibration/` evaluates limited risk-output updating in hospitals held out
from the updating fit. These stages must not be pooled or renamed as one
analysis.

The public R scripts are path-adapted copies of the as-run scripts. They use
`EHR_AUDIT_REPO_ROOT`, `EHR_AUDIT_RESTRICTED_DATA_ROOT`, and `EHR_AUDIT_OUTPUT_ROOT`.
The restricted root contains `final_full.csv`, `baseline_covars.csv`, and
`eicu_external.csv` directly. The output root contains `model/` with the
locked partition objects and the archived `domain3_deployable_generator_v1_20260728/`
analysis, and receives the new derived objects and private predictions. This
is an explicit local contract, not a claim that the scripts run on an arbitrary
PhysioNet download layout.
`EHR_AUDIT_OUTPUT_ROOT` must exist outside the repository; each script checks this
before writing restricted objects. The output paths are ignored by Git, but
ignoring them is not a substitute for keeping them outside the checkout.
Run recovery, L1 assignment, L1 model fitting, L1 bootstrap, and finally
hospital-held-out recalibration in that order. The original hash-checked
recalibration plan remains in the controlled workspace and is supplied by
`D3_RECAL_ORIGINAL_PLAN`; the public copy does not invent an equivalent
pre-run hash. The recovered parameter object is not distributed pending a
separate data-governance decision. This is therefore not a credential-free
end-to-end replay.

Files in a `private_local` output subdirectory are restricted intermediates
and must never be committed. This release contains no D3 empirical output CSVs.
# Current source version: 2026-10-08

`source_correction/` contains the current post hoc temperature-corrected L1
external path. Other directories retain earlier prerequisites and historical
comparators; their external outputs are not automatically corrected.
Only Base + K2 external predictions are updated by correction script 03.
The subsequent `raw_temperature_completion/` restores and verifies source
fits, retains the archived L2 Raw recipe, evaluates corrected external Raw
inputs for both outcomes and supplies current paired comparisons/Figure 4C.
The historical fitting scripts did not persist their original Raw fits;
the completion outputs verified reconstructions, not retrieved original fits.
See `docs/D3_CALIBRATION_SOURCE_NOTE.md` and the source crosswalk.

