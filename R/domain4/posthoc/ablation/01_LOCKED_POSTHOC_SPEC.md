# D4 feature-block ablation comparison: locked post hoc specification

Written on 2026-09-25 (Asia/Shanghai), after the original K2 D4 findings but before generating the ablation labels or computing any new D4 comparison result. The author's approval of one new common, partition-agnostic propensity model supersedes the earlier instruction to reuse K2-conditioned patient scores and weights. The original D4 model and results remain untouched.

## Scientific question and claim limit

Compare measured treatment-support diagnostics under the locked K2 partition and one Step 0 feature-block ablation partition, using the same 17,525-patient retrospective hour-24 risk set, the same first-RRT-positive-derived-record indicator, one common propensity model, and the same patient-level weights. The only variable changed **within the new comparison** is group membership. This post hoc comparison cannot establish that removing a feature block causally changes treatment assignment, that a record is an active RRT start, or that an RRT effect or effect modification is identified.

## Fixed source and label rule

- Source full cohort: `final_full.csv`; fixed K2: `labels_primary_mice.rds`; standardized completed matrix: `X_primary33_std_mice.rds`; formal D4 risk-set ledger: `D4_landmark_patient_ledger_private.rds`. Join exclusively by unique `stay_id`; do not rely on row order after imports.
- Step 0 ablation is **A4: remove the registry's strict renal-core plus strict acid-base blocks** from the 33-feature standardized matrix, then run Lloyd K-means, `K=2`, `nstart=25`, `iter.max=100`, `algorithm="Lloyd"`. Its 20-seed output saved summary metrics but no located patient-level label vector. Therefore choose the originally specified first seed, `20240601`. Reproduce the archived ARI for this seed before proceeding; do not substitute another seed. Also reproduce locked K2 using the all-33-feature A0 first-seed check.
- Verify that the A4 input features are exactly the 33-variable list excluding the seven features registered as strict renal/acid-base: `urine_output_24h_ml`, `bun_max`, `creatinine_max`, `ph_min`, `pco2_max`, `aniongap_max`, and `bicarbonate_min`. Record every retained feature. Other related variables, including lactate, potassium and chloride, remain; describe the object accordingly.
- Orient the two A4 labels using the archived `final_full$mortality_30d` on the entire 20,049-patient cohort: higher observed death proportion is A1, lower A2. If tied, stop rather than choose an arbitrary orientation. This is a label-name operation, not a basis for selecting a seed or optimizing D4 support.
- Report full-cohort and D4 risk-set sizes and the full 2x2 K2-by-A4 cross-tabulation. Stop if either A4 group has fewer than 500 patients in the risk set.

## One common propensity model

- Keep the exact original risk set, outcome and event-order fields. Preserve archived K2 score, ATE weight and published support tables as separate historical results. Do not use their K2-conditioned weights in the new between-partition comparison.
- Fit exactly one new pooled, unpenalized binomial `glm` for the derived first-RRT-positive-record indicator. Covariates are the original D4 list **except** `cluster_k2`: age, gender, first-day SOFA, first-day peak AKI stage, SAPS II, maximum creatinine, maximum BUN, maximum potassium, minimum bicarbonate, minimum pH, maximum lactate, first-day urine total, minimum mean blood pressure, maximum heart rate, maximum respiratory rate, and minimum GCS. No A4 label enters the model.
- Reuse the original D4 median/mode imputation algorithm, `gender` factor handling, logistic link, and probability clipping to `[1e-6, 1-1e-6]`. Record every imputation fill, warnings, convergence, coefficient finiteness, probability range, clipping count and model formula. Stop on nonconvergence or nonfinite coefficients. Do not tune the covariates, link, clip or patient subset after seeing diagnostics.
- Compute for each patient from this one score `e`: original-scale ATE weights `A/e+(1-A)/(1-e)` and overlap weights `A(1-e)+(1-A)e`. Use identical patient-level values under K2 and A4. Retain scores and weights only in `private_not_for_release` and never publish or upload patient-level rows.

## Original-result checks before new comparison

- Assert risk-set `N=17,525`, locked K2 counts `3,151/14,374`, first-record-positive counts `265/49`, and archived death-cell counts `149/992` and `18/1,762`. The first-record source and 30-day mortality producer limitations documented in the earlier D4 package remain unresolved.
- Check that the new model has exactly the expected 16 non-partition covariates and no partition label. Its weights are new; the old K2 ESS and balance verdict are **not** numerical reproduction targets. Separately verify that the archived K2 propensity file remains unchanged. For the new K2 and A4 diagnostic, compute both from the same new score and the same balance code.

## Metrics and decision screen

For K2 and A4, each group separately: risk-set count, positive-record count/rate, positive and negative cell archived 30-day death/missing counts; ATE positive/negative-arm ESS, ESS/arm n, top-ceiling-1%-of-positive-arm normalized-weight share; empirical min/max score-range limits and positive/negative arm outside counts; all available prespecified covariate standardized differences before weighting, after ATE and after overlap weighting. Standardize every weighted difference by the group's fixed **unweighted pooled arm SD**. Record zero/undefined denominators rather than applying an epsilon. The four negative source urine totals make urine balance data-quality affected: display it, but exclude it from the decision screen. Report the full 15 interpretable variable list and every absolute SMD above 0.1, especially kidney and acid-base variables. Compute between-group unweighted SMDs for the seven strict block variables plus retained lactate/potassium/chloride as a mechanism-premise description.

A group meets the **post hoc measured-balance/event screen** if every one of the 15 interpretable covariates has overlap-weighted absolute SMD **strictly below 0.1** on the unrounded scale and each record-positive/negative cell has at least 10 observed archived 30-day deaths. The original 1% positive-record frequency threshold is reported but not included in this new screen. Both groups meeting the screen is the between-group screen result. Do not label this causal identification, clinical equipoise, or evidence of RRT treatment heterogeneity.

Prediction fixed now: both A4 groups meet the screen. If not, retain the failure and report exact group(s) and covariate(s), particularly kidney/acid-base balance. No alternative seed, group boundary, propensity model, covariate deletion, SMD denominator or threshold may be selected in response.

## Output and privacy

Write a source/selection report, aggregate 2x2 cross-tabs, a long aggregate metrics table, a threshold summary, and an English Methods/Results draft of at most 120 words. The manuscript and supplement DOCX remain unchanged. Save a private label/score audit object only under `private_not_for_release`; it is not a release artifact. No RRT effect, interaction, mortality regression, new outcome horizon, or patient-level public file is authorized.

## Source fingerprints at lock

| Source | SHA-256 |
|---|---|
| `X_primary33_std_mice.rds` | `9896AD4906207CA910E54EB187507037106BFF212299586141908F387CFE6279` |
| `labels_primary_mice.rds` | `159576D0B2796E9786DB0EB4C8747AF987BE5E53AF200419E3E37567C1D875F9` |
| `final_full.csv` | `CE0B38D2D3D6263A5C4615A6AEA1E78B19840BD4B35046DBFA65EE6BFBBA633B` |
| `D4_landmark_patient_ledger_private.rds` | `233EBC8593B6789B04A7F84C66BDBC19BBCEC503D08573D4892C8FAB7C0486A7` |
| `38_d4_landmark_support_common_v1.R` | `2B6ECC405E15924CB69A77C7ACB3A2E88DE79D2E8DE99D397B0C9759DF1F78DD` |
| `Table_D1_ablation_seed_stability.csv` | `888DE81D3C67897D0A501907938B249F060B2D731830D09548B12CCCFB02B6FC` |
| `D1_feature_domain_registry.csv` | `A20D2EAC3C698654E5724E3AE5C1C49BAABC58C1C08E72C401963D8894631673` |
