-- Bounded aggregate comparison of archived eICU GCS with current raw rows.
BEGIN TRANSACTION READ ONLY;
SET LOCAL statement_timeout = '180s';

WITH raw_gcs AS MATERIALIZED (
  SELECT n.patientunitstayid, n.nursingchartvalue::numeric AS value_numeric
  FROM nursecharting n
  WHERE n.nursingchartoffset BETWEEN 0 AND 1440
    AND lower(coalesce(n.nursingchartcelltypevallabel, ''))
        LIKE '%glasgow coma score%'
    AND lower(coalesce(n.nursingchartcelltypevalname, '')) LIKE '%gcs total%'
    AND n.nursingchartvalue ~ '^[0-9]+$'
), source_gcs AS (
  SELECT r.patientunitstayid, min(r.value_numeric) AS gcs_min
  FROM raw_gcs r
  JOIN public.eicu_external e USING (patientunitstayid)
  GROUP BY r.patientunitstayid
)
SELECT count(*) AS cohort_stays,
       count(e.gcs_min) AS archived_gcs_observed,
       count(g.gcs_min) AS current_raw_gcs_observed,
       count(*) FILTER (WHERE e.gcs_min IS NULL AND g.gcs_min IS NOT NULL)
         AS current_raw_only,
       count(*) FILTER (WHERE e.gcs_min IS NOT NULL AND g.gcs_min IS NULL)
         AS archived_only,
       count(*) FILTER (WHERE e.gcs_min IS NOT NULL AND g.gcs_min IS NOT NULL
                         AND e.gcs_min <> g.gcs_min) AS both_observed_disagree,
       percentile_cont(0.5) WITHIN GROUP (ORDER BY abs(e.gcs_min - g.gcs_min))
         AS median_abs_gcs_difference
FROM public.eicu_external e
LEFT JOIN source_gcs g USING (patientunitstayid);

ROLLBACK;
