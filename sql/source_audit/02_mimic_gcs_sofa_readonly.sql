-- Aggregate-only MIMIC source/timing audit. Run in the authorized MIMIC-IV
-- database. The first_day_gcs definition must also be checked against the
-- exact local mimic-code version; a table value alone cannot reveal it.
BEGIN TRANSACTION READ ONLY;

SELECT table_schema, table_name, column_name
FROM information_schema.columns
WHERE table_schema = 'mimiciv_derived'
  AND table_name IN ('first_day_gcs', 'sofa', '_metadata')
  AND column_name IN ('gcs_min', 'gcs_unable', 'gcs_verbal',
                      'starttime', 'endtime', 'sofa_24hours',
                      'attribute', 'value')
ORDER BY table_name, ordinal_position;

SELECT count(*) AS cohort_stays,
       count(f.gcs_min) AS observed_first_day_gcs,
       count(*) FILTER (WHERE f.gcs_min = 15) AS gcs_15_stays,
       count(*) FILTER (WHERE f.gcs_min BETWEEN 3 AND 14) AS gcs_3_to_14_stays
FROM project_sa_aki.strict_sepsis3_cohort c
LEFT JOIN mimiciv_derived.first_day_gcs f USING (stay_id);

SELECT count(*) AS cohort_stays,
       count(sofa_time) AS assessed_stays,
       count(*) FILTER (WHERE sofa_time > intime + interval '24 hours')
         AS sofa_assessed_after_icu_hour_24,
       count(*) FILTER (WHERE sofa_time < intime) AS sofa_assessed_before_icu
FROM project_sa_aki.strict_sepsis3_cohort;

ROLLBACK;
