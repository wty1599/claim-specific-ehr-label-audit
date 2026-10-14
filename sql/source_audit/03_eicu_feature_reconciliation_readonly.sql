-- Compare archived eICU features with raw first-day nurseCharting values.
-- Aggregate output only; no cohort rows are exported.
BEGIN TRANSACTION READ ONLY;

WITH nursing AS (
  SELECT n.patientunitstayid,
         n.nursingchartoffset,
         lower(coalesce(n.nursingchartcelltypevallabel, '')) AS label,
         lower(coalesce(n.nursingchartcelltypevalname, '')) AS name,
         CASE WHEN btrim(coalesce(n.nursingchartvalue, '')) ~
                '^[+-]?[0-9]+(\.[0-9]+)?$'
              THEN btrim(n.nursingchartvalue)::numeric END AS value_numeric
  FROM nursecharting n
  JOIN public.eicu_external e USING (patientunitstayid)
  WHERE n.nursingchartoffset BETWEEN 0 AND 1440
    AND lower(coalesce(n.nursingchartcelltypevallabel, '')) = 'temperature'
    AND lower(coalesce(n.nursingchartcelltypevalname, '')) IN
        ('temperature (c)', 'temperature (f)')
), valid AS (
  SELECT patientunitstayid, nursingchartoffset, name,
         CASE WHEN name = 'temperature (c)' AND value_numeric BETWEEN 25 AND 45
                   THEN value_numeric
              WHEN name = 'temperature (f)' AND value_numeric BETWEEN 77 AND 113
                   THEN (value_numeric - 32) * 5 / 9 END AS celsius
  FROM nursing
), per_stay AS (
  SELECT patientunitstayid,
         min(celsius) FILTER (WHERE name = 'temperature (c)') AS nurse_c_min,
         min(celsius) FILTER (WHERE name = 'temperature (f)') AS nurse_f_min,
         min(celsius) AS nurse_any_min,
         count(celsius) AS valid_rows
  FROM valid
  GROUP BY patientunitstayid
)
SELECT count(*) AS cohort_stays,
       count(e.temperature_min) AS periodic_observed,
       count(p.nurse_c_min) AS nurse_c_observed,
       count(p.nurse_f_min) AS nurse_f_observed,
       count(p.nurse_any_min) AS nurse_any_observed,
       count(*) FILTER (WHERE e.temperature_min IS NULL AND p.nurse_any_min IS NOT NULL)
         AS newly_observed_by_nursing,
       count(*) FILTER (WHERE e.temperature_min IS NOT NULL AND p.nurse_any_min IS NOT NULL)
         AS both_sources,
       percentile_cont(0.5) WITHIN GROUP
         (ORDER BY abs(e.temperature_min - p.nurse_any_min)) AS median_abs_min_difference_c,
       count(*) FILTER (WHERE abs(e.temperature_min - p.nurse_any_min) > 1)
         AS min_difference_over_1_c,
       count(*) FILTER (WHERE abs(e.temperature_min - p.nurse_any_min) > 2)
         AS min_difference_over_2_c
FROM public.eicu_external e
LEFT JOIN per_stay p USING (patientunitstayid);

WITH nursing AS (
  SELECT n.patientunitstayid, n.nursingchartoffset,
         lower(coalesce(n.nursingchartcelltypevalname, '')) AS name,
         CASE WHEN btrim(coalesce(n.nursingchartvalue, '')) ~
                   '^[+-]?[0-9]+(\.[0-9]+)?$'
              THEN btrim(n.nursingchartvalue)::numeric END AS value_numeric
  FROM nursecharting n
  JOIN public.eicu_external e USING (patientunitstayid)
  WHERE n.nursingchartoffset BETWEEN 0 AND 1440
    AND lower(coalesce(n.nursingchartcelltypevallabel, '')) = 'temperature'
    AND lower(coalesce(n.nursingchartcelltypevalname, '')) IN
        ('temperature (c)', 'temperature (f)')
), paired AS (
  SELECT patientunitstayid, nursingchartoffset,
         avg(value_numeric) FILTER (WHERE name = 'temperature (c)'
             AND value_numeric BETWEEN 25 AND 45) AS c,
         avg((value_numeric - 32) * 5 / 9) FILTER
             (WHERE name = 'temperature (f)' AND value_numeric BETWEEN 77 AND 113) AS f_to_c
  FROM nursing
  GROUP BY patientunitstayid, nursingchartoffset
)
SELECT count(*) FILTER (WHERE c IS NOT NULL AND f_to_c IS NOT NULL)
         AS same_offset_c_f_pairs,
       count(DISTINCT patientunitstayid) FILTER
         (WHERE c IS NOT NULL AND f_to_c IS NOT NULL) AS stays_with_pairs,
       percentile_cont(0.5) WITHIN GROUP (ORDER BY abs(c - f_to_c))
         AS median_abs_unit_discrepancy_c,
       count(*) FILTER (WHERE abs(c - f_to_c) > 0.5) AS pairs_over_half_c
FROM paired;

SET LOCAL statement_timeout = '180s';

WITH raw_gcs AS MATERIALIZED (
  SELECT n.patientunitstayid,
         n.nursingchartvalue::numeric AS value_numeric
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
