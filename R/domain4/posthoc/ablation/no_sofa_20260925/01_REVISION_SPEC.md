# D4 ablation comparison: post hoc no-SOFA sensitivity specification

Written 2026-09-25 (Asia/Shanghai) after the original D4 result **and after the first partition-neutral comparison**. The first comparison and all earlier source objects remain unchanged. This second version is authorized by the author in response to a source-timing audit, not prespecified before the original study and not selected to improve a support verdict. Hash this document before running any revised model.

## Question and interpretation boundary

On the unchanged retrospective 17,525-patient hour-24 risk set, compare K2 and the single fixed A4 feature-block-ablation partition under **one new common logistic model of the first RRT-positive record before competing events**. Remove `sofa_score`, whose qualification timestamp was not retained in `final_full` and for which a parallel, otherwise exactly matching export placed 151 records after hour 24. The first version, with SOFA, is a required side-by-side sensitivity comparator. Removing SOFA from a model does not remove possible post-landmark *cohort eligibility*: the matched 151 records suggest some Sepsis-3 qualification occurred later. Moreover, SAPS II uses first-day physiology but its standard MIMIC implementation includes comorbidities inferred from hospitalization diagnosis codes; it cannot be asserted to be fully available at hour 24. The local joined SAPS II total has not been decomposed, and this version retains it because the author authorized removal of SOFA only. Other first-day source tables use chart time, with an inclusive hour-24 boundary, rather than a verified result-availability time; the KDIGO urine component can inherit retrospectively backfilled weight. These are unresolved temporal provenance limitations, not proof that every retained value was available when a clinician could act. Neither version identifies an RRT treatment effect, treatment heterogeneity, clinical equipoise, or a prospective hour-24 decision cohort.

## Allowed change, protected objects and privacy

- Allowed change: remove **only** `sofa_score` from the v1 partition-neutral propensity covariates, from the weighted-balance variable list and from the resulting measured-balance screen. The new model thus has 15 original non-partition covariates, and the decision screen has 14 interpretable variables after also excluding the source-quality-affected urine total. Do not delete or change any other variable.
- Protected: original K2 labels and all D1-D3 results; A4 first seed 20240601, label orientation, full and risk-set membership; hour-24 RRT-free/surviving/in-ICU risk set; 24-to-72-hour first-event order and same 30-day mortality; v1 formula/results; original D4 formula/results; no model fitting for mortality or RRT effects.
- Fit exactly one new pooled, unpenalized binomial GLM, with age, gender, first-day AKI stage, SAPS II, maximum creatinine/BUN/potassium/lactate, minimum bicarbonate/pH/mean BP/GCS, first-day urine total, maximum heart rate and respiratory rate. No K2 or A4 term, no SOFA term, and no mortality outcome term. Reuse median/mode imputation, gender factor handling, logit link, probability clipping `[1e-6, 1-1e-6]`, and all unchanged RRT record definitions. The inherited SAPS II coding caveat remains. Record formula, imputation fills, warnings, convergence, raw/clipped score range and coefficients' finiteness. Stop for nonconvergence or nonfinite coefficients.
- Compute patient-level ATE and overlap weights once from the single new score, then use exactly the same values under K2 and A4. Store all patient-level labels, scores and weights only under this version's `private_not_for_release`; never place them in a public package.

## Fixed source checks

Before fitting, verify source hashes, unique stay IDs and 1:1 joins. Assert risk-set N=17,525; K2 C1/C2 n=3,151/14,374; RRT-positive counts=265/49; archived 30-day death cells C1 positive/negative=149/992 and C2=18/1,762; A4 group counts=5,168/12,357; A4 vs K2 seed/ARI and label source hash match v1. Verify the v1 aggregate and private objects remain unchanged before and after the run. Audit the remaining 15 covariates' documented source windows and explicitly distinguish chart time from information-availability time. Record SAPS II diagnosis-code, first-day inclusive-boundary, urine-output/weight, and unresolved `final_full` assembly caveats; do not silently remove a second covariate. If source inspection proves that another retained field directly includes observations strictly after hour 24, stop and report rather than revise the model again. Unverified hour-24 availability is a claim limitation for this narrowly authorized sensitivity run.

## Measures and display

For each partition group report: n, first RRT-positive record n/rate, positive and negative cells' observed 30-day deaths/missingness, ATE/overlap arm-specific ESS and ESS/n, maximum normalized weight and top-ceiling-1%-of-positive-arm weight share, empirical score-range overlap limits and arm-specific outside counts, and before/after-weighting SMDs for all 15 retained model covariates. Every weighted SMD uses that partition group's **fixed unweighted pooled-arm SD**, unchanged across weight types. Do not add epsilons to undefined SDs or structural zero cells. Report all variables exceeding absolute 0.1, especially kidney/acid-base terms. The first version's SOFA row is not eligible for the fair 14-variable comparison.

Retain the previously written **post hoc descriptive screen**, without claiming causal support: a group meets it only if every one of 14 interpretable (non-urine) retained variables has overlap-weighted absolute SMD strictly `<0.1` at full precision **and** both record-positive and record-negative cells have at least 10 observed archived 30-day deaths. Both groups must meet it for a partition-level descriptive screen result. Report the original >=1% record-positive frequency criterion separately; it is not in the screen. No new thresholds, partitions, covariates, seeds, probability truncations or treatment models may be chosen after viewing results.

Compare v1 and v2 on the **same 14 variables**, along with arm-specific ESS and record-positive counts. V1's original 15-variable screen can be shown separately but not mislabelled as a same-variable contrast. The original A4-both-groups-pass prediction, fixed before v1, is retained as a historical prediction; v2 is an outcome-aware source-timing correction and therefore cannot be called its independent confirmatory test. Report regardless of whether balance improves or worsens.

## Claim ceiling and deliverables

The model estimates which patients had a first recorded RRT-positive event before death, ICU exit or administrative end, not treatment assignment at hour 24. Competing events and variable observation opportunity remain in the modelled record outcome. Even if an overlap-weighted subgroup has low measured SMDs, it is not a causal RRT effect or proof of positivity. Both partitions share one retrospective risk set; they define different subgroups, not different population risk sets.

Write only new-version scripts, private scores, aggregate results, source manifests, QA and a Chinese report with an English manuscript-ready **conditional** draft. Do not edit the submission manuscript or supplement in this task. If the temporal/provenance audit does not pass, deliver the failure report instead of a clean manuscript paragraph.

## Anchored sources from v1

| Source | Locked SHA-256 |
|---|---|
| Parent original D4 risk ledger | `233EBC8593B6789B04A7F84C66BDBC19BBCEC503D08573D4892C8FAB7C0486A7` |
| Parent A4 patient-level label vector | `0A0E9AFAB6D320B48CBD5C7A13B814A46C5C3F8A84DF47DC3A69763A13795987` |
| Original D4 imputation audit | `D85600EC7BFB2E4D6C949BD1E2FF432C0B47FDF76091885349A5FC89B0A3A019` |
| Parent v1 patient-level common-score object | `CA2B096AE452776BFD243F2116A2036A64F92ACB9BC5BEDCA72BA7B70E459887` |
| Parent v1 balance CSV | `923C6D2511F92750DC69CF3B62E2C72337AB40F13AC37C7B0FB4D2BC10893B59` |
| Parent v1 group support CSV | `164113E5F7AE7EA617E6B64D04399A445A4F2FF03CA417EFDFE5049C4E809954` |
| Parent v1 group screen CSV | `2A770AAE55F0A7363BC75F1FC50184267FD0188FEB4C27AB639D3D84C94D1602` |
| Parent v1 corrected between-group screen CSV | `1D8B55AFBB810523A66C9E33522C8D4556B1790DFF2EABFE31560AC763498343` |
| Parent v1 post-landmark SOFA aggregate | `84F59CD655FB47A26BF005C05E2C179D992AA6040D03394BB5AB57769BBDA14A` |
