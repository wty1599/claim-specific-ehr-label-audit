-- Post hoc external-temperature source export under the frozen specification.
-- Read-only database transaction; per-stay output stays in private_local.
BEGIN TRANSACTION READ ONLY;
SET LOCAL statement_timeout = '300s';
\pset format csv
-- Supply temperature_export_path as a psql variable pointing outside this repository.
\o :temperature_export_path

WITH nursing AS MATERIALIZED (
  SELECT n.patientunitstayid,
         lower(coalesce(n.nursingchartcelltypevalname, '')) AS name,
         CASE WHEN btrim(coalesce(n.nursingchartvalue, '')) ~
                     '^[+-]?[0-9]+(\.[0-9]+)?$'
              THEN btrim(n.nursingchartvalue)::numeric END AS value_numeric
  FROM nursecharting n
  WHERE n.nursingchartoffset BETWEEN 0 AND 1440
    AND lower(coalesce(n.nursingchartcelltypevallabel, '')) = 'temperature'
    AND lower(coalesce(n.nursingchartcelltypevalname, '')) IN
        ('temperature (c)', 'temperature (f)')
), nursing_stay AS (
  SELECT n.patientunitstayid,
         min(n.value_numeric) FILTER
           (WHERE n.name = 'temperature (c)' AND n.value_numeric BETWEEN 25 AND 45)
             AS nurse_c_min,
         max(n.value_numeric) FILTER
           (WHERE n.name = 'temperature (c)' AND n.value_numeric BETWEEN 25 AND 45)
             AS nurse_c_max,
         min((n.value_numeric - 32) * 5 / 9) FILTER
           (WHERE n.name = 'temperature (f)' AND n.value_numeric BETWEEN 77 AND 113)
             AS nurse_f_min_c,
         max((n.value_numeric - 32) * 5 / 9) FILTER
           (WHERE n.name = 'temperature (f)' AND n.value_numeric BETWEEN 77 AND 113)
             AS nurse_f_max_c
  FROM nursing n
  JOIN public.eicu_external e USING (patientunitstayid)
  GROUP BY n.patientunitstayid
)
SELECT e.patientunitstayid,
       CASE WHEN n.nurse_c_min IS NOT NULL THEN n.nurse_c_min
            WHEN n.nurse_f_min_c IS NOT NULL THEN n.nurse_f_min_c
            WHEN e.temperature_min BETWEEN 25 AND 45 THEN e.temperature_min END
         AS temperature_min_corrected,
       CASE WHEN n.nurse_c_max IS NOT NULL THEN n.nurse_c_max
            WHEN n.nurse_f_max_c IS NOT NULL THEN n.nurse_f_max_c
            WHEN e.temperature_max BETWEEN 25 AND 45 THEN e.temperature_max END
         AS temperature_max_corrected,
       CASE WHEN n.nurse_c_min IS NOT NULL THEN 'nurse_c'
            WHEN n.nurse_f_min_c IS NOT NULL THEN 'nurse_f_converted'
            WHEN e.temperature_min BETWEEN 25 AND 45
              OR e.temperature_max BETWEEN 25 AND 45 THEN 'periodic_fallback'
            ELSE 'missing' END AS temperature_source
FROM public.eicu_external e
LEFT JOIN nursing_stay n USING (patientunitstayid)
ORDER BY e.patientunitstayid;

\o
ROLLBACK;
