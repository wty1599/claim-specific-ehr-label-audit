# SQL extraction guide

The SQL was run against local PostgreSQL installations of MIMIC-IV and eICU.
Database access, credentialing, and the MIMIC-IV derived tables are not provided
by this repository. Execute with `psql` and stop on errors.

## MIMIC-IV order

1. `mimic/00_build_adult_first_icu_cohort.sql`
2. `mimic/01_strict_sepsis3.sql`
3. `mimic/02_kdigo_aki_outcome.sql`
4. `mimic/03_outcomes.sql`
5. `mimic/04_export_final_sepsis_analysis.sql`
6. `mimic/05_extract_features_24h.sql`
7. `mimic/06_build_baseline_covars.sql`
8. `mimic/07_assemble_initial_final_full.sql`
9. `mimic/08_study_specific_renal_composite.sql`
10. `mimic/09_merge_renal_composite.sql`

Scripts `00` and `07` are recovered production blocks. The original standalone
files were absent from the audited current source roots; their SQL logic was
recovered from the prior submission-code release and requires author
verification before a public release is tagged.

The final submitted MIMIC export had 20,049 rows and 86 columns. Script `08`
creates a study-specific 30-day renal composite. It must not be described as a
canonical MAKE30 endpoint.

Export paths in psql scripts can be overridden with `-v`, for example:

```text
psql ... -v final_full_csv='<restricted-data-root>/final_full.csv' -f sql/mimic/09_merge_renal_composite.sql
```

## eICU order

1. `eicu/01_build_proxy_sepsis_cohort_mdrd.sql`
2. `eicu/02_extract_features_24h.sql`
3. `eicu/03_build_modified_sofa.sql`
4. `eicu/04_deduplicate_and_export.sql`
5. `eicu/05_flowchart_counts.sql` for descriptive counts

The eICU sepsis definition uses diagnosis-string matching and is a proxy rather
than a direct implementation of the MIMIC Sepsis-3 onset algorithm. The
submitted export contained 17,465 admissions from 199 hospitals after prior
CKD/ESRD exclusions. The physical table name `public.saaki_cohort` is historical
and retained to avoid breaking downstream SQL.

## Domain 4 extraction

Run `domain4/01_extract_rrt_timing_inputs.sql` after the final MIMIC table is
available. The four public output-path variables are:

- `rrt_first_csv`
- `rrt_windows_csv`
- `rrt_base_times_csv`
- `features_0_6_csv`

These exports are patient-level and must remain outside the public repository.

## Encoding note

Distributed SQL files use English comments, UTF-8 without a byte-order mark,
and LF line endings. Comment translation did not change executable SQL tokens.

