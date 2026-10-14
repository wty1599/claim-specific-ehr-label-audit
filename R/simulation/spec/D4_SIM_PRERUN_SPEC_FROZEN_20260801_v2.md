# Frozen specification v2: longitudinal D4 support-interface validation

Frozen on 2026-08-01 before the v2 smoke and formal simulations. Version 1 is
retained unchanged as an audit artifact. Version 2 clarifies the competing-event
generator and adds source-hash, patient-level equivalence, and threshold-boundary
requirements. Audit thresholds and expected control states are unchanged.

## Aim and estimand boundary

Validate the current D4 support-diagnostic interface under known longitudinal
data-generating mechanisms. The target is correct classification of whether the
data support model development. No treatment effect is generated, estimated, or
interpreted.

## Fixed audit interface

- Time zero: ICU admission; landmark: ICU hour 24.
- Opportunity risk set: alive, in the index ICU, and RRT-free before hour 24.
- Exposure: RRT as the first event during [24, 72) hours.
- Competing/precluding events: death and ICU exit.
- Exact-tie priority: RRT, death, ICU exit, administrative hour 72.
- Label: fixed K2 partition generated before the landmark; no refitting.
- Propensity model: fixed K2 label plus the locked 16 empirical covariates.
- Missing covariates: the empirical median/mode interface.
- Support metrics are calculated separately within C1 and C2: outside common
  support, untruncated ATE-IPW ESS/N, minimum treatment cell, and proportion of
  weights greater than 10.

## Locked gates and thresholds

- Primary Gate B: treated proportion >=0.01 and >=10 30-day mortality events in
  both treated and untreated cells within each phenotype.
- Outside-overlap alert: proportion >0.10.
- ESS/N alert: ratio <0.25.
- Minimum-cell alert: count <30.
- Extreme-weight alert: proportion with ATE weight >10 exceeds 0.01.
- Any Gate B failure or phenotype-specific support alert yields
  `SUPPORT_INADEQUATE`; otherwise the state is
  `SUPPORT_ADEQUATE_FOR_MODEL_DEVELOPMENT`.

Threshold equality is non-alerting. No threshold, truth label, coefficient, or
scenario definition may be changed after source hashes are frozen.

## Shared covariate and label mechanism

For each patient, severity S~N(0,1), renal R=0.55S+sqrt(1-0.55^2)e_R, and shock
H=0.50S+sqrt(1-0.50^2)e_H. The fixed partition score is
0.85R+0.55H+0.25e_K; C1 is the upper 20% and C2 the remainder. The empirical D4
covariates are deterministic noisy functions of S, R, and H as implemented in
`39b_d4_sim_generate_longitudinal_controls_v1.R`. EHR-like MAR missingness is
applied only after event generation.

## Pre-landmark mechanism

In G4-L and G6-L, pre-ICU RRT probability is 0.0002; pre-24-hour RRT, death, and
ICU exit probabilities use the fixed logistic formulas in the generator. D4-LQ
uses probabilities 0.010, 0.080, 0.080, and 0.180, respectively, to stress the
exclusion ledger.

## Independent candidate-time generator

Eligible patients receive independent RRT, death, and ICU-exit candidate times.
For a cause probability p over 48 hours, the exponential rate is
`-log(1-p)/48`; waits >=48 hours are administratively censored. The ledger, not
the generator, selects the first event. This prevents sequential cause sampling
from mechanically determining event order.

### G4-L: longitudinal support-failure positive control

The RRT cause probability uses phenotype-specific baselines of 0.10 (C1) and
0.0045 (C2), with coefficients 0.65 for within-phenotype centered renal status
and 0.25 for within-phenotype centered shock. Expected state:
`SUPPORT_INADEQUATE`.

### G6-L: longitudinal adequate-support negative control

The RRT cause probability is `logit^-1(logit(0.22)+0.18S+0.12R+0.10I[C1])`,
bounded to [0.18, 0.35]. Expected state:
`SUPPORT_ADEQUATE_FOR_MODEL_DEVELOPMENT`.

### D4-LQ: ledger stress control

D4-LQ uses the G6-L allocation mechanism bounded to [0.12, 0.40], increases
pre-landmark exclusions, and imposes exact RRT/death ties in approximately 1%
of eligible patients. It is excluded from operating-characteristic denominators.

For all scenarios, death cause probability is
`logit^-1(-3.75+0.55S+0.35I[C1])`; ICU-exit cause probability is
`logit^-1(-0.55-0.35S)`. A generated death before hour 72 always implies
30-day mortality; additional deaths through day 30 follow the fixed outcome
model in the generator.

## Simulation size

- Formal initial cohort: 20,049; 300 repetitions each for G4-L and G6-L; 100
  repetitions for D4-LQ.
- Smoke cohort: 2,000; 5 repetitions each for G4-L and G6-L; 3 for D4-LQ.
- Base random seed: 20260801.

## Mandatory qualification checks

1. All executable source hashes must match the separately frozen v2 manifest.
2. Reapplication to the locked empirical source must exactly reproduce the
   hour-24 patient set, first-event states, event times, tie flags, treatment
   indicator, event counts, support metrics, Gate A-D table, and claim state.
3. Unit tests must confirm inclusive/exclusive behavior at every locked threshold.
4. Failures must be attributed to generation, ledger, propensity/imputation, or
   metric/state stages and reported under exclusion and failure-as-error rules.
5. All D4-LQ repeats must contain exact ties resolved using the locked priority.

## Performance measures

Sensitivity for G4-L, specificity for G6-L, scenario-specific correct
classification, Wilson 95% confidence intervals, continuous-metric mean, SD,
empirical 2.5th-97.5th percentiles and Monte Carlo SE, and failure audits.

## Prohibited analyses

No treatment effect, HTE, interaction, risk difference, risk ratio, odds ratio,
E-value, or causal estimand is fitted or reported.
