# Completion of corrected-temperature external Raw model evaluation

Recorded: 2026-10-08. This is a post hoc completion authorized after the
temperature-source correction and its Base+K2 results were known.

The original Raw-EN and Raw-RF fit objects were not retained. Reconstruct the
original fits using the archived development inputs, preprocessing recipe,
model families, tuning procedure, iteration order, random seed, and runtime.
Do not change the MIMIC partition, outcomes, covariates, cohort, hyperparameters,
or external outcome-evaluable subsets. Do not run clustering or MICE.

Before evaluating corrected inputs, require the reconstructed models to
reproduce every archived external prediction on the original extract, aligned
by outcome, model, and patientunitstayid. Tolerance: 1e-12. Stop on failure;
do not relabel a new fit as the original fit.

Apply the verified fits and their unchanged original Raw-model preprocessing
to the corrected eICU temperature minimum and maximum. All other external
inputs remain byte-identical in memory. Use the existing mortality (17,284
stays, 2,940 events) and descriptive composite (17,333 stays, 5,714 events)
subsets. Report AUC, Brier, joint calibration slope/intercept and patient-level
percentile bootstrap AUC intervals, 1,000 replicates. The current Base+K2
predictions and hospital-held-out recalibration are not refitted.

Complete current paired external AUC contrasts against Base+K2 using shared
patient bootstrap indices, 1,000 replicates. Record the seed, failures and
conditional fixed-prediction uncertainty. Preserve historical results as such;
write new results in an independent directory. Restricted predictions and
recovered fit objects remain local. Only code and documentation can enter the
new public repository.

Synchronize the deterministic Figure 4C data, rendering, alternative text and
affected current external-performance supplement rows. Do not change Figure
4A/B/D data or D1, D2, D4 results. A substantial change in the main D3 conclusion
requires author review before narrative revision.
