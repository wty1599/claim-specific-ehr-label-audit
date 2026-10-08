-- Aggregate-only audit. Run in the authorized eICU database with the source
-- tables on search_path and public.eicu_external from the locked analysis.
-- Do not export patient-level rows or publish rare source-label cells.
BEGIN TRANSACTION READ ONLY;

SELECT count(*) AS cohort_stays,
       count(temperature_min) AS periodic_temperature_stays,
       count(gcs_min) AS nurse_total_gcs_stays
FROM public.eicu_external;

WITH rows AS (
  SELECT n.patientunitstayid,
         lower(coalesce(n.nursingchartcelltypecat, '')) AS category,
         lower(coalesce(n.nursingchartcelltypevallabel, '')) AS label,
         lower(coalesce(n.nursingchartcelltypevalname, '')) AS name,
         btrim(coalesce(n.nursingchartvalue, '')) AS raw_value
  FROM nursecharting n
  JOIN public.eicu_external c USING (patientunitstayid)
  WHERE n.nursingchartoffset BETWEEN 0 AND 1440
), temperature AS (
  SELECT *,
         CASE WHEN raw_value ~ '^[+-]?[0-9]+(\.[0-9]+)?$'
              THEN raw_value::numeric END AS value_numeric
  FROM rows
  WHERE label LIKE '%temp%' OR name LIKE '%temp%'
)
SELECT category, label, name,
       count(*) AS rows,
       count(DISTINCT patientunitstayid) AS stays,
       count(value_numeric) AS numeric_rows,
       count(*) FILTER (WHERE value_numeric BETWEEN 25 AND 45) AS celsius_range_rows,
       count(*) FILTER (WHERE value_numeric BETWEEN 77 AND 113) AS fahrenheit_range_rows
FROM temperature
GROUP BY category, label, name
ORDER BY stays DESC, category, label, name;

WITH temperature AS (
  SELECT DISTINCT n.patientunitstayid
  FROM nursecharting n
  JOIN public.eicu_external c USING (patientunitstayid)
  WHERE n.nursingchartoffset BETWEEN 0 AND 1440
    AND (lower(coalesce(n.nursingchartcelltypevallabel, '')) LIKE '%temp%'
      OR lower(coalesce(n.nursingchartcelltypevalname, '')) LIKE '%temp%')
    AND btrim(coalesce(n.nursingchartvalue, '')) ~ '^[+-]?[0-9]+(\.[0-9]+)?$'
)
SELECT count(*) AS cohort_stays,
       count(c.temperature_min) AS periodic_temperature_stays,
       count(t.patientunitstayid) AS any_numeric_nurse_temperature_stays,
       count(*) FILTER (WHERE c.temperature_min IS NULL AND t.patientunitstayid IS NOT NULL)
         AS nurse_only_candidate_stays,
       count(*) FILTER (WHERE c.temperature_min IS NOT NULL AND t.patientunitstayid IS NOT NULL)
         AS both_source_stays
FROM public.eicu_external c
LEFT JOIN temperature t USING (patientunitstayid);

WITH rows AS (
  SELECT n.patientunitstayid,
         lower(coalesce(n.nursingchartcelltypevallabel, '')) AS label,
         lower(coalesce(n.nursingchartcelltypevalname, '')) AS name,
         btrim(coalesce(n.nursingchartvalue, '')) AS raw_value
  FROM nursecharting n
  JOIN public.eicu_external c USING (patientunitstayid)
  WHERE n.nursingchartoffset BETWEEN 0 AND 1440
    AND lower(coalesce(n.nursingchartcelltypevallabel, '')) LIKE '%glasgow coma score%'
)
SELECT label, name, count(*) AS rows,
       count(DISTINCT patientunitstayid) AS stays,
       count(*) FILTER (WHERE raw_value ~ '^[0-9]+$') AS integer_rows,
       count(*) FILTER (WHERE CASE WHEN raw_value ~ '^[0-9]+$'
          THEN raw_value::numeric END BETWEEN 3 AND 15) AS valid_total_rows
FROM rows
GROUP BY label, name
ORDER BY stays DESC, label, name;

SELECT count(*) AS cohort_stays,
       count(c.gcs_min) AS nurse_total_gcs_stays,
       count(DISTINCT a.patientunitstayid) FILTER
         (WHERE a.eyes IS NOT NULL AND a.motor IS NOT NULL AND a.verbal IS NOT NULL)
         AS complete_apache_components,
       count(*) FILTER
         (WHERE c.gcs_min IS NOT NULL AND a.eyes IS NOT NULL
           AND a.motor IS NOT NULL AND a.verbal IS NOT NULL)
         AS both_gcs_sources
FROM public.eicu_external c
LEFT JOIN apacheapsvar a USING (patientunitstayid);

ROLLBACK;
