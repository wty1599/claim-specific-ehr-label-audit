# Post hoc D3 external temperature source correction

Fixed before computing any revised external labels or outcome metrics.
Date: 2026-10-08 (Asia/Shanghai). This is a response to the source audit
after the initial D3 results were known, not part of the original protocol.

## Inputs and protected objects

- Population: the same 17,465 eICU stays in the archived L1 extract.
- Keep the 20,049 fixed MIMIC-IV K=2 labels, recovered 33-feature recipe,
  scaling parameters, and centroids unchanged. Do not recluster.
- Keep the archived eICU extract and all previous D3 outputs read-only.
- Keep all non-temperature external features, including the existing
  `gcs_min` and modified SOFA, unchanged in this correction. GCS measurement
  differences between databases are disclosed separately; APACHE components
  are not substituted for charted GCS totals.
- No outcome, original prediction, or hospital-level performance is used to
  choose the temperature rule.

## Temperature rule

Use ICU minutes 0 through 1440, inclusive, as in the archived feature SQL.
For each stay, accept numeric `nursecharting` rows with category label
`Temperature` and value name `Temperature (C)` in 25-45 degrees C. The
minimum and maximum of those values are the primary first-day temperature
features. If no valid C row exists, accept `Temperature (F)` in 77-113
degrees F, convert each value to C using `(F-32)*5/9`, and take its minimum
and maximum. If neither nursing source exists, retain a finite
`vitalperiodic` minimum/maximum only when it lies in 25-45 degrees C.
Otherwise the feature is missing and receives the original fixed-source
median through the unchanged L1 recipe. Do not use `Temperature Location`
as a temperature measurement.

The nursing C and F entries usually describe the same event. They are not
added as independent observations. For stays with both nursing and periodic
values, nursing is primary because MIMIC first-day temperature was drawn
from charted vitals; the source difference is documented rather than
silently averaged. Record overlap, missingness, and discordance in aggregate.

## Analyses and reporting

1. Produce one new restricted 17,465-row temperature source file with the
   unchanged stay identifiers and corrected min/max. Assert unique stays,
   source coverage, valid ranges, and equality of all non-temperature
   feature values with the archived extract.
2. Apply the unchanged recovered L1 rule to the modified external extract.
   Report label changes, direction, prevalence, missingness strata, and
   hospital aggregate variation. Do not inspect outcomes before this step.
3. If labels change, rerun the original D3 mortality and non-equivalent
   composite models with the same covariates, splits, fitting and interval
   procedures, followed by the same hospital-held-out recalibration and
   disclosure-aware Figure 4 source data. Do not tune to recover old numbers.
4. If any main D3 conclusion changes, report the revised source-audit and
   aggregate model results before changing the manuscript. Keep old and new
   results visibly distinct. The manuscript and public code are not edited
   by the source export itself.

All patient-level labels, feature rows, and predictions remain in the
authorized local `private_local` area and are not uploaded.
