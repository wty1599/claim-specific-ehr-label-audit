# Post hoc D3 hospital-prevalence analysis specification

Fixed on 2026-10-04 (Asia/Shanghai), after the original D3 prevalence and
release-disclosure results were known and before the funnel analysis below.
This is not an original prespecified analysis. The original 69/199 statistic
motivated the analysis. Do not select the center, limits, display rule, or
interpretation by looking at the new results.

Publication adaptation dated 2026-10-08: private source and document paths
are replaced by input roles below. This preserves the original 2026-10-04
statistical plan. For the corrected run, the temperature-correction
specification supersedes the historical external counts: scripts 06 and 07
use corrected_c1, total C1=3319, and the exact external center 3319/17465.
The original reference center, limits, model and display rules remain fixed.

## Scope and immutable inputs

- Analyze only one row per eICU hospital: L1 assigned-stay count `n_stays`
  and L1 C1 count `l1_c1_n`. No individual records or patient-level data are
  read. Hospital identifiers may be used only inside the restricted analysis
  for the random intercept and are never exported.
- Historical restricted input role: hospital counts for the uncorrected L1 rule.
- Corrected restricted input role: hospital counts produced by script 02,
  stored outside this repository. Its fixed input SHA-256 is checked in
  scripts 06 and 07. Neither table is distributed.
- Primary fixed reference: MIMIC-IV L1 C1 prevalence `3902/20049`.
  Sensitivity reference: eICU L1 C1 prevalence `3313/17465`. Display these
  as 19.46% and 18.97%, but use the exact fractions for calculations.

## Source and disclosure checks

Require 199 unique hospitals, integer counts with `0 <= C1 <= n`, total
`n=17465`, total `C1=3313`, and 69 hospitals whose observed prevalence
differs by more than 0.10 from `3902/20049`. A failure stops the run.

The release-review rule is **not a single total-n threshold**: a hospital is
eligible for an individual public point only if it has at least 10 C1 stays
and at least 10 C2 stays. The previous review yielded 96 individually shown
and 103 suppressed hospitals. Verify these counts. The distinction matters
because a large hospital with a small C1 cell may still be suppressed.

## Calculations

1. Report the distribution of `n_stays` across 199 hospitals: minimum,
   quartiles, median, maximum, selected size bands, and the 96/103 split.
   Split the previously reported 69 hospitals with >10 percentage-point
   difference by individual-display eligibility.
2. For each hospital n and each fixed center p, form exact central equal-tail
   binomial control limits with integer counts
   `L=qbinom(alpha/2,n,p)` and `U=qbinom(1-alpha/2,n,p)`.
   A hospital is below if `C1<L`, above if `C1>U`, otherwise within.
   Calculate separately for alpha 0.05 (95%) and 0.002 (99.8%).
   Report above, below, and total out of 199, both overall and by display
   eligibility. Discreteness makes actual tail probabilities no larger than
   the nominal limits. These are unadjusted descriptive limits, not a
   multiple-testing correction or case-mix-adjusted test.
3. Fit the single planned grouped-binomial logistic random-intercept model
   `cbind(C1,n-C1) ~ 1 + (1|hospital)` using `lme4::glmer`, `nAGQ=9`, and
   `bobyqa` with `maxfun=200000`. Report the hospital random-intercept SD on
   the logit scale and its 95% profile-likelihood CI. Check convergence and
   singularity; do not swap optimizers or fit a second model to obtain a
   desired result. This model accounts for binomial sampling under its
   assumptions, but the variance can still reflect case mix and repeated
   stays, so it is not a pure biological hospital effect.

## Interpretation stop

Treat zero, one, or two hospitals beyond the MIMIC-centered 99.8% limits as
"few" for the purpose of the requested editorial stop. This was fixed before
calculation because 199 x 0.002 is about 0.4 nominal false alarms; discrete
limits are often more conservative. If the count is at most two, report the
results and ask the author to choose revised hospital-language wording before
editing the manuscript or figures. Also stop if the random-intercept model
does not converge or is singular, if its profile CI cannot be obtained, or if
any source/disclosure assertion fails. Counts above two permit figure and
limited text revision but do not by themselves establish causal or case-mix-
independent hospital heterogeneity.

## Figure 4B and public source if the stop is not triggered

- Plot hospital n on a log x-axis and L1 C1 prevalence on y. Show the
  MIMIC-centered 95% and 99.8% binomial limits. Draw and directly label
  the MIMIC and eICU pooled prevalence lines with different line types.
- Public individual points are limited to hospitals with both C1 and C2
  counts at least 10. Strip hospital IDs and avoid stable rank keys.
- Aggregate the other hospitals in predetermined n bands: 1-4, 5-9,
  10-19, 20-39, 40-79, and 80 or more stays. Merge adjacent bands from low
  to high, if needed, until each released group contains at least two
  hospitals, at least 10 total stays, and at least 10 C1 and 10 C2 stays.
  If the last band fails, merge backward. Group prevalence is total C1/total
  stays; its plotted x coordinate is the median hospital n and is explicitly
  identified as a grouped summary, not an individual hospital. If safe bands
  cannot be formed, stop. Do not export subgroup hospital IDs or member rows.
- No output is to be added to a code repository in this task. The 199-row
  restricted source remains private. Any source CSV released with the figure
  must be reviewed against this disclosure rule.

## Reporting scope

Distinguish 17,465 assignable stays from 17,284 mortality-evaluable and
17,333 composite-evaluable stays. The funnel analysis does not establish
semantic equivalence, treatment prediction, or a causal hospital effect.
