# D4 landmark support prerun specification

Status: **FROZEN - AUTHORIZED FOR FORMAL EXECUTION**

Frozen: 2026-07-30  
Authorization basis: author instruction to implement and archive the D4 plan.

## Locked scientific object

- Cohort: 20,049 MIMIC-IV stays.
- Labels: `output/model/labels_primary_mice.rds`.
- Expected K2 counts: C1=3,992; C2=16,057.
- No clustering or relabeling.

## Locked time interface

- Landmark: ICU hour 24.
- Eligible at landmark:
  - `outtime - intime >= 24`;
  - death missing or `death_time_hours >= 24`;
  - first RRT missing or `first_rrt_hours >= 24`.
- Treatment window: first RRT in `[24,72)` hours.
- Observation end: earliest of RRT, death, ICU exit, or hour 72.
- Primary exact-tie priority: RRT, death, ICU exit, hour 72.
- Mandatory tie sensitivity: exclude tied rows.

## Locked event states

- `RRT_FIRST_24_72`
- `DEATH_FIRST_24_72`
- `ICU_EXIT_FIRST_24_72`
- `EVENT_FREE_AT_72`

## Locked propensity formula

`rrt_first_24_72 ~ cluster_k2 + age + gender + sofa_score +
aki_stage_0_24h + sapsii + creatinine_max + bun_max + potassium_max +
bicarbonate_min + ph_min + lactate_max + urine_output_24h_ml + mbp_min +
heart_rate_max + resp_rate_max + gcs_min`

- Primary missing handling: median/mode imputation, with values exported.
- Primary estimator: unpenalized binomial logistic regression.
- No fallback estimator may replace the primary silently.
- No outcome or post-hour-24 covariate is allowed.

## Locked treatment-opportunity gates

- A: treated fraction >=0.005 and mortality events >=5 in treated and untreated.
- B: treated fraction >=0.010 and mortality events >=10 in treated and untreated.
- C: treated fraction >=0.020 and mortality events >=10 in treated and untreated.
- D: treated fraction >=0.010 and mortality events >=20 in treated and untreated.
- Primary gate: B.
- Mortality field for the event-count screen: `mortality_30d` from
  `data/final_full.csv`.

## Locked structural support alerts

Calculated separately in C1 and C2:

- outside-overlap proportion >0.10;
- ATE-IPW ESS/N <0.25;
- minimum treated/untreated cell <30;
- proportion with ATE-IPW >10 exceeds 0.01.

## Locked state machine

1. Missing interface or failed structural QC -> `INTERFACE_INCOMPLETE`.
2. Complete interface plus any primary Gate B failure in either K2 group ->
   `SUPPORT_INADEQUATE`.
3. Complete interface plus any structural support alert in either K2 group ->
   `SUPPORT_INADEQUATE`.
4. Complete interface and no Gate B or structural support failure ->
   `SUPPORT_ADEQUATE_FOR_MODEL_DEVELOPMENT`.

## Prohibited outputs

- treatment effect;
- HTE or treatment-by-K2 interaction;
- weighted RD, RR, or OR;
- E-value;
- causal-benefit language.

## Formal run authorization

- Scientific specification: FROZEN
- Input hashes: FROZEN in the execution plan and common configuration
- Formal output: new versioned directory only
- Overwrite of WP10 or submission outputs: PROHIBITED

