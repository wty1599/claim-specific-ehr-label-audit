# Post hoc RRT support diagnostics

`01_compute_empirical_support.R` reads the original hour-24 risk-set ledger
and saved propensity audit from the restricted workspace, checks the original
17,525-patient identity, and emits aggregate-only arm-specific ESS, weight
concentration, overlap, fixed-scale balance, and target-population summaries.
It does not refit the propensity model or estimate an RRT effect.

The public copy requires `D4_FORMAL_OUTPUT_ROOT` to point to the archived
`output_formal_v1` directory and uses `EHR_AUDIT_OUTPUT_ROOT` for new results. The latter receives
local outputs under `domain4_arm_support/empirical_v2`. The source RDS files
contain patient records and are not distributed. Reviewed copies of the
aggregate outputs are under `data/empirical_results/domain4_arm_support/`.
`EHR_AUDIT_OUTPUT_ROOT` must exist outside the repository; the script checks this
before writing. The published ESS table omits individual weight extrema and
small-cell top-weight diagnostics retained only in restricted output.
