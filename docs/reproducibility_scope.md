# Reproducibility Scope

## Local Preparation On 2026-10-08

Version `v1.0.0-jamia-submission` (2026-10-08) contains previously public Git-tracked sources and specifically authorized correction scripts, SQL and dated specifications, without inherited Git history or outputs. Publication adaptations affect paths, commentary, messages and private authorization prerequisites, not numerical algorithms.

No database query, clinical analysis, model fit, bootstrap campaign, missing-pattern transfer or figure render was run during preparation. No restricted row-level input was read. Source/destination hashes and review records remain outside the code directory. Static parsing and artificial-fixture tests are not numerical reproduction of clinical results.

## Correction Scope

The 2026-10-08 specification is post hoc: written after earlier D3 findings were known and before its own calculations. It corrects first-day temperature for the same 17,465 assignable eICU stays and applies the unchanged recovered L1 rule. MIMIC-IV labels and centroids are not re-estimated.

Script 03 recalculates **only Base + K2 external predictions** with archived fitted coefficients. Script 04 bootstraps those predictions and retains archived comparators. Script 05 repeats mortality recalibration with corrected Base + K2 and archived Base, preserving patient-linked hospital components and held-out allocation. Scripts 06-08 prepare the corrected funnel and Figure 4 sources. Their supplied code was not executed during this release preparation.

Base contains no temperature feature. Raw-EN and Raw-RF use Base covariates plus 33 inputs including temperature. The subsequent authorized Raw completion reconstructs source fits because original fit objects were not saved, verifies all four historical models for both outcomes within 1e-12, and evaluates the two Raw models on corrected external temperature with the unchanged archived L2 recipe. It does not substitute the recovered L1 assignment scale. Current Figure 4C and paired comparisons use completed corrected Raw sources; correction 03/04 remain the earlier Base + K2-only steps, not the final current Raw comparison.

Model recovery/evaluation requires the original R 4.6.0 and asserted package versions; rendering retains R 4.5.3. The source execution record reports successful old-prediction agreement and zero failures in 1,000 patient-bootstrap replicates per current statistic. Packaging did not repeat model fitting, bootstrap, patient-level comparisons or rendering. The separate completion verifier is a deterministic source/summary check, not an independent model rerun.

Completion 03 preserves fully reproduced historical estimates/intervals when CSV rounding perturbs RF rank ties. The calibration-definition addendum retains offset-only CITL and the joint-fit slope, with the joint intercept separately named. Step 04 stages current panel C and unchanged A/B/D inputs by byte-preserving copy; no plot calculation or empirical value is recomputed.

Script 11 encodes the separately frozen post hoc transfer of corrected eICU missingness masks to 4,689 complete MIMIC-IV records: 500 repetitions, seed 20261008, no reclustering or outcome model. This estimates technical assignment movement, not external accuracy, semantic validity or the causal fraction of a prevalence difference. It was not run here.

## Historical Limits

Clinical extracts, matrices, labels, predictions, fitted models, parameter objects, hospital rows, aggregate outputs and figure datasets are excluded. Exact replay requires authorized access and retained intermediate inputs. Hash assertions for restricted scientific inputs are retained; public specification hashes identify distributed specification bytes, not private prompts or authorization records.

Original D2 ten-fold assignments were not saved. Public code cannot restore them. Historical D1 controls, complete simulation campaigns, the original cohort-wide MICE workflow and source SQL were not rerun here. Excluded historical registries and result objects remain execution prerequisites, not evidence of an independent validation.

The clinical illustration is one sepsis partition. It does not empirically establish wider transportability of the specification. Hospital prevalence and calibration do not establish class discreteness, causal hospital effects, external semantic equivalence or treatment prediction. The loop-diuretic feasibility check stopped at its timing screen; no treatment-effect models followed.

## Publication Adaptations

`EHR_AUDIT_*` configuration and relative code paths replace workstation locations. English strings and consistent historical folder aliases replace Chinese strings. Distributed public counterparts replace available code dependencies. Private work-package authorization and prompt checks are excluded from the public contract; removal details remain in the separate local audit. Unavailable scientific inputs are not fabricated.

Main-figure input and output paths are external. The shared figure adapter gives an explicit Figure 4 argument precedence over its environment/default root. Renderers record the files they actually consume; the entry verifies those records against unchanged selected-source fingerprints. These path checks were tested with artificial files only.

One current Figure 4 renderer preserves the authoritative final plotting logic and assertions. Required preparation scripts remain, redundant older preparation versions are removed, and empirical sources/pixels are excluded. Local code preparation was limited to source packaging and static/artificial-fixture checks.
