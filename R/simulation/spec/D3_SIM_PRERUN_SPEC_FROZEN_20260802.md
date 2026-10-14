# Frozen specification: regenerated-generator D3 validation

Frozen before the formal run. Smoke testing may correct executable defects but
must not change the data-generating mechanisms, truth labels, thresholds, or
formal sample sizes below. Any executable change after smoke requires a new
source-hash manifest before formal execution.

## Aim

Validate the current D3 computational lineage under known data-generating
mechanisms: a MICE-imputation-1 Lloyd K=2 discovery partition, a separately
regenerated deterministic source-median nearest-centroid generator, frozen
external application, and a Baseline-plus-generated-K2 mortality model.

The simulation does not claim that the regenerated generator is the exact
original discovery generator. Internal and external fidelity are reported
against explicit reference objects.

## Alignment with the locked semi-synthetic harness

- 33 EHR-like first-24-hour physiological features.
- Formal derivation n=5,000 and external n=3,000.
- 100 repetitions per scenario.
- MICE m=5 and maxit=5; completed dataset 1 defines the locked partition.
- Formal 1st/99th percentile winsorization after MICE, derivation-set
  standardization, and Lloyd K=2 k-means with 50 starts.
- C1 is oriented as the higher-severity group, as in the legacy simulation.
- EHR-like MAR missingness is 8% in derivation and non-drift external data.
- The legacy G3 covariate-drift external cohort uses severity mean 0.55,
  severity SD 1.10, the six frozen feature shifts, and MNAR-like missingness
  intensity 18%.

## Regenerated generator

The generator is built from the incomplete derivation data and the locked K=2
labels only. For each feature it freezes observed-data 1st/99th percentiles,
the post-winsor median, and observed post-winsor mean/SD. Centroids are means
within the locked labels in this surrogate coordinate system. New records are
assigned by squared Euclidean nearest centroid. No external outcome, external
preprocessing refit, or external reclustering is used.

## Scenarios and truth

1. `G3R_same_distribution_no_alert`: derivation and external cohorts share the
   same feature, missingness, and outcome mechanisms. Combined D3 alert truth
   is negative.
2. `G3R_covariate_drift_alert`: the external cohort reproduces the legacy G3
   severity, feature, and missingness drift. Combined D3 alert truth is
   positive.
3. `D3CQ_calibration_stress_control`: features and missingness remain
   same-distribution, but the external mortality generator changes its
   intercept from -2.35 to -1.85 and its severity coefficient from 0.90 to
   1.35. This is a calibration-engine positive control. Combined D3 alert truth
   is positive. It is reported separately from the legacy covariate-drift
   control.

## Locked D3 rules

- Pooled prevalence alert: absolute C1 prevalence drift >0.10.
- Assignment-fidelity alert: regenerated assignment versus the simulated
  complete-feature frozen reference ARI <0.80.
- Calibration alert: Baseline-plus-generated-K2 mortality calibration slope
  outside [0.80, 1.20] or absolute calibration intercept >0.20.
- Combined D3 alert: any component alert.
- Equality at a threshold is non-alerting.

Internal generator agreement and ARI are continuous fidelity results. No
retrospective pass threshold is introduced for internal fidelity.

## Reference objects

- Internal reference: the MICE-imputation-1 Lloyd K=2 locked label vector.
- External simulation-only reference: the discovery preprocessing and
  centroids applied to the complete pre-missingness external features. This is
  available only because the data-generating process is known; it is not an
  empirical eICU reference partition.

## Formal and smoke sizes

- Formal: n=5,000 derivation, n=3,000 external, 100 repetitions per scenario.
- Smoke: n=1,000 derivation, n=600 external, 2 repetitions per scenario.
- Base seed: 20260802.
- Formal worker default: 4 on the 16-GB host; workers may be reduced without
  changing the random-number registry.

## Failure handling and summaries

Generation, MICE, discovery, generator construction, assignment, and outcome
model failures are recorded by stage. The primary operating-characteristic
summary excludes failed repeats and reports the failure rate. A second summary
counts failures as incorrect classifications. Proportions use Wilson 95%
intervals. Continuous metrics use mean, SD, empirical 2.5th/97.5th percentiles,
and Monte Carlo SE=SD/sqrt(n).

## Prohibited changes

- No threshold or truth label may be tuned after inspecting results.
- No external outcome may enter generator construction or assignment.
- No external reclustering may determine the generated label.
- No old G3 output is overwritten or renamed as a new result.
- No manuscript or supplementary document is modified in this run.
