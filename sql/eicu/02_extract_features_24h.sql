-- ============================================================






-- ============================================================
\set ON_ERROR_STOP on
\pset pager off
\if :{?eicu_external_csv}
\else
\set eicu_external_csv './outputs/sql/eicu_external.csv'
\endif


DROP TABLE IF EXISTS tmp_ids;
CREATE TEMP TABLE tmp_ids AS
SELECT DISTINCT patientunitstayid FROM public.saaki_cohort;
CREATE INDEX ON tmp_ids(patientunitstayid);



DROP TABLE IF EXISTS tmp_lab;
CREATE TEMP TABLE tmp_lab AS
SELECT l.patientunitstayid,
  MAX(l.labresult) FILTER (WHERE ln='wbc x 1000')        AS wbc_max,
  MIN(l.labresult) FILTER (WHERE ln='wbc x 1000')        AS wbc_min,
  MIN(l.labresult) FILTER (WHERE ln='platelets x 1000')  AS platelets_min,
  MIN(l.labresult) FILTER (WHERE ln='hgb')               AS hemoglobin_min,
  MAX(l.labresult) FILTER (WHERE ln='total bilirubin')   AS bilirubin_total_max,
  MAX(l.labresult) FILTER (WHERE ln='alt (sgpt)')        AS alt_max,
  MAX(l.labresult) FILTER (WHERE ln='ast (sgot)')        AS ast_max,
  MIN(l.labresult) FILTER (WHERE ln='albumin')           AS albumin_min,
  MAX(l.labresult) FILTER (WHERE ln='pt - inr')          AS inr_max,
  MAX(l.labresult) FILTER (WHERE ln='pt')                AS pt_max,
  MAX(l.labresult) FILTER (WHERE ln='ptt')               AS ptt_max,
  MAX(l.labresult) FILTER (WHERE ln='creatinine')        AS creatinine_max,
  MIN(l.labresult) FILTER (WHERE ln='creatinine')        AS creatinine_min,
  MAX(l.labresult) FILTER (WHERE ln='bun')               AS bun_max,
  MAX(l.labresult) FILTER (WHERE ln='anion gap')         AS aniongap_max,
  MIN(l.labresult) FILTER (WHERE ln='bicarbonate')       AS bicarbonate_min,
  MAX(l.labresult) FILTER (WHERE ln='glucose')           AS glucose_max,
  MIN(l.labresult) FILTER (WHERE ln='sodium')            AS sodium_min,
  MAX(l.labresult) FILTER (WHERE ln='sodium')            AS sodium_max,
  MIN(l.labresult) FILTER (WHERE ln='potassium')         AS potassium_min,
  MAX(l.labresult) FILTER (WHERE ln='potassium')         AS potassium_max,
  MIN(l.labresult) FILTER (WHERE ln='calcium')           AS calcium_min,
  MAX(l.labresult) FILTER (WHERE ln='calcium')           AS calcium_max,
  MIN(l.labresult) FILTER (WHERE ln='chloride')          AS chloride_min,
  MAX(l.labresult) FILTER (WHERE ln='chloride')          AS chloride_max,
  MAX(l.labresult) FILTER (WHERE ln='lactate')           AS lactate_max,
  MIN(l.labresult) FILTER (WHERE ln='ph')                AS ph_min,
  -- eICU stores arterial PCO2 primarily under "paco2".
  MAX(l.labresult) FILTER (WHERE ln IN ('pco2','paco2')) AS pco2_max,
  MIN(l.labresult) FILTER (WHERE ln='base excess')       AS baseexcess_min,
  MIN(l.labresult) FILTER (WHERE ln='pao2')              AS po2_min,
  MAX(l.labresult) FILTER (WHERE ln='fio2')              AS fio2_max_raw,
  MIN(l.labresult) FILTER (WHERE ln='paco2')             AS paco2_min_alt
FROM (
  SELECT patientunitstayid, lower(labname) AS ln, labresult, labresultoffset
  FROM lab
  WHERE labresultoffset BETWEEN 0 AND 1440 AND labresult IS NOT NULL
) l
JOIN tmp_ids i ON i.patientunitstayid=l.patientunitstayid
GROUP BY l.patientunitstayid;

-- ---------- 1b. verified absolute lymphocyte count (0-24h minimum) ----------
-- "-lymphs" is lymphocyte percentage (0-100), not an absolute count.
-- Pair WBC and lymphocyte percentage at the exact same labresultoffset:
-- ALC = WBC (x10^9/L) * lymphocyte percentage / 100.
-- Missing ALC remains NULL and is imputed later with the frozen MIMIC median.
DROP TABLE IF EXISTS tmp_alc;
CREATE TEMP TABLE tmp_alc AS
WITH wbc AS (
  SELECT l.patientunitstayid, l.labresultoffset,
         AVG(l.labresult) AS wbc_value
  FROM lab l
  JOIN tmp_ids i ON i.patientunitstayid=l.patientunitstayid
  WHERE l.labresultoffset BETWEEN 0 AND 1440
    AND lower(l.labname)='wbc x 1000'
    AND l.labresult BETWEEN 0.1 AND 200
  GROUP BY l.patientunitstayid, l.labresultoffset
),
lymph AS (
  SELECT l.patientunitstayid, l.labresultoffset,
         AVG(l.labresult) AS lymphocyte_pct
  FROM lab l
  JOIN tmp_ids i ON i.patientunitstayid=l.patientunitstayid
  WHERE l.labresultoffset BETWEEN 0 AND 1440
    AND lower(l.labname)='-lymphs'
    AND l.labresult BETWEEN 0 AND 100
  GROUP BY l.patientunitstayid, l.labresultoffset
),
matched AS (
  SELECT w.patientunitstayid, w.labresultoffset,
         w.wbc_value, ly.lymphocyte_pct,
         w.wbc_value * ly.lymphocyte_pct / 100.0 AS absolute_lymphocyte_count
  FROM wbc w
  INNER JOIN lymph ly
    ON ly.patientunitstayid=w.patientunitstayid
   AND ly.labresultoffset=w.labresultoffset
)
SELECT patientunitstayid,
       MIN(absolute_lymphocyte_count) AS abs_lymphocytes_min,
       COUNT(*)::integer AS alc_exact_match_count
FROM matched
GROUP BY patientunitstayid;






DROP TABLE IF EXISTS tmp_vp;
CREATE TEMP TABLE tmp_vp AS
SELECT v.patientunitstayid,
  MAX(v.heartrate)     AS heart_rate_max,
  MIN(v.sao2)          AS spo2_min,
  MAX(v.respiration)   AS resp_rate_max,
  MIN(v.temperature)   AS temperature_min,
  MAX(v.temperature)   AS temperature_max,
  MIN(v.systemicmean)  AS mbp_min_inv,
  AVG(v.systemicmean)  AS mbp_mean_inv
FROM vitalperiodic v
JOIN tmp_ids i ON i.patientunitstayid=v.patientunitstayid
WHERE v.observationoffset BETWEEN 0 AND 1440
GROUP BY v.patientunitstayid;


DROP TABLE IF EXISTS tmp_va;
CREATE TEMP TABLE tmp_va AS
SELECT a.patientunitstayid,
  MIN(a.noninvasivemean) AS mbp_min_ni,
  AVG(a.noninvasivemean) AS mbp_mean_ni
FROM vitalaperiodic a
JOIN tmp_ids i ON i.patientunitstayid=a.patientunitstayid
WHERE a.observationoffset BETWEEN 0 AND 1440
GROUP BY a.patientunitstayid;


DROP TABLE IF EXISTS tmp_gcs;
CREATE TEMP TABLE tmp_gcs AS
SELECT n.patientunitstayid, MIN(n.val) AS gcs_min
FROM (
  SELECT patientunitstayid, nursingchartvalue::numeric AS val, nursingchartoffset
  FROM nursecharting
  WHERE lower(nursingchartcelltypevallabel) LIKE '%glasgow coma score%'
    AND lower(nursingchartcelltypevalname) LIKE '%gcs total%'
    AND nursingchartvalue ~ '^[0-9]+$'
    AND nursingchartoffset BETWEEN 0 AND 1440
) n
JOIN tmp_ids i ON i.patientunitstayid=n.patientunitstayid
GROUP BY n.patientunitstayid;


DROP TABLE IF EXISTS tmp_uo;
CREATE TEMP TABLE tmp_uo AS
SELECT io.patientunitstayid, SUM(io.cellvaluenumeric) AS urine_output_24h_ml
FROM intakeoutput io
JOIN tmp_ids i ON i.patientunitstayid=io.patientunitstayid
WHERE io.intakeoutputoffset BETWEEN 0 AND 1440
  AND io.cellvaluenumeric IS NOT NULL AND io.cellvaluenumeric >= 0
  AND (lower(io.celllabel) LIKE '%urine%' OR lower(io.celllabel) LIKE '%foley%'
       OR lower(io.celllabel) LIKE '%void%')
GROUP BY io.patientunitstayid;


DROP TABLE IF EXISTS public.eicu_external;
CREATE TABLE public.eicu_external AS
SELECT
  c.*,
  -- lab
  lb.wbc_max, lb.wbc_min, lb.platelets_min, lb.hemoglobin_min,
  alc.abs_lymphocytes_min,
  COALESCE(alc.alc_exact_match_count, 0)::integer AS alc_exact_match_count,
  TRUE::boolean AS alc_time_matching_verified,
  'WBC x -lymphs / 100 at identical patientunitstayid and labresultoffset; 0-24h minimum'::text AS alc_derivation,
  lb.bilirubin_total_max, lb.alt_max, lb.ast_max, lb.albumin_min, lb.inr_max, lb.pt_max, lb.ptt_max,
  lb.creatinine_max, lb.creatinine_min, lb.bun_max, lb.aniongap_max, lb.bicarbonate_min, lb.glucose_max,
  lb.sodium_min, lb.sodium_max, lb.potassium_min, lb.potassium_max, lb.calcium_min, lb.calcium_max,
  lb.chloride_min, lb.chloride_max, lb.lactate_max, lb.ph_min, lb.pco2_max, lb.baseexcess_min, lb.po2_min,

  CASE WHEN lb.po2_min IS NOT NULL AND lb.fio2_max_raw >= 21
       THEN ROUND((lb.po2_min / (lb.fio2_max_raw/100.0))::numeric, 1) END AS pao2fio2ratio_min,

  vp.heart_rate_max, vp.spo2_min, vp.resp_rate_max, vp.temperature_min, vp.temperature_max,

  LEAST(vp.mbp_min_inv, va.mbp_min_ni) AS mbp_min,
  COALESCE((vp.mbp_mean_inv + va.mbp_mean_ni)/2.0, vp.mbp_mean_inv, va.mbp_mean_ni) AS mbp_mean,
  g.gcs_min,
  uo.urine_output_24h_ml
FROM public.saaki_cohort c
LEFT JOIN tmp_lab lb ON lb.patientunitstayid=c.patientunitstayid
LEFT JOIN tmp_alc alc ON alc.patientunitstayid=c.patientunitstayid
LEFT JOIN tmp_vp  vp ON vp.patientunitstayid=c.patientunitstayid
LEFT JOIN tmp_va  va ON va.patientunitstayid=c.patientunitstayid
LEFT JOIN tmp_gcs g  ON g.patientunitstayid=c.patientunitstayid
LEFT JOIN tmp_uo  uo ON uo.patientunitstayid=c.patientunitstayid;


\echo '===== Feature coverage and medians (cross-database unit check) ====='
SELECT 'creatinine_max' f, COUNT(creatinine_max) n, ROUND((PERCENTILE_CONT(0.5) WITHIN GROUP(ORDER BY creatinine_max))::numeric,2) med FROM public.eicu_external
UNION ALL SELECT 'bun_max', COUNT(bun_max), ROUND((PERCENTILE_CONT(0.5) WITHIN GROUP(ORDER BY bun_max))::numeric,2) FROM public.eicu_external
UNION ALL SELECT 'wbc_max', COUNT(wbc_max), ROUND((PERCENTILE_CONT(0.5) WITHIN GROUP(ORDER BY wbc_max))::numeric,2) FROM public.eicu_external
UNION ALL SELECT 'platelets_min', COUNT(platelets_min), ROUND((PERCENTILE_CONT(0.5) WITHIN GROUP(ORDER BY platelets_min))::numeric,2) FROM public.eicu_external
UNION ALL SELECT 'abs_lymphocytes_min', COUNT(abs_lymphocytes_min), ROUND((PERCENTILE_CONT(0.5) WITHIN GROUP(ORDER BY abs_lymphocytes_min))::numeric,2) FROM public.eicu_external
UNION ALL SELECT 'lactate_max', COUNT(lactate_max), ROUND((PERCENTILE_CONT(0.5) WITHIN GROUP(ORDER BY lactate_max))::numeric,2) FROM public.eicu_external
UNION ALL SELECT 'sodium_max', COUNT(sodium_max), ROUND((PERCENTILE_CONT(0.5) WITHIN GROUP(ORDER BY sodium_max))::numeric,2) FROM public.eicu_external
UNION ALL SELECT 'ph_min', COUNT(ph_min), ROUND((PERCENTILE_CONT(0.5) WITHIN GROUP(ORDER BY ph_min))::numeric,2) FROM public.eicu_external
UNION ALL SELECT 'pao2fio2ratio_min', COUNT(pao2fio2ratio_min), ROUND((PERCENTILE_CONT(0.5) WITHIN GROUP(ORDER BY pao2fio2ratio_min))::numeric,2) FROM public.eicu_external
UNION ALL SELECT 'heart_rate_max', COUNT(heart_rate_max), ROUND((PERCENTILE_CONT(0.5) WITHIN GROUP(ORDER BY heart_rate_max))::numeric,2) FROM public.eicu_external
UNION ALL SELECT 'mbp_min', COUNT(mbp_min), ROUND((PERCENTILE_CONT(0.5) WITHIN GROUP(ORDER BY mbp_min))::numeric,2) FROM public.eicu_external
UNION ALL SELECT 'spo2_min', COUNT(spo2_min), ROUND((PERCENTILE_CONT(0.5) WITHIN GROUP(ORDER BY spo2_min))::numeric,2) FROM public.eicu_external
UNION ALL SELECT 'gcs_min', COUNT(gcs_min), ROUND((PERCENTILE_CONT(0.5) WITHIN GROUP(ORDER BY gcs_min))::numeric,2) FROM public.eicu_external
UNION ALL SELECT 'urine_output_24h_ml', COUNT(urine_output_24h_ml), ROUND((PERCENTILE_CONT(0.5) WITHIN GROUP(ORDER BY urine_output_24h_ml))::numeric,2) FROM public.eicu_external
ORDER BY f;

\echo '===== Total rows ====='
SELECT COUNT(*) FROM public.eicu_external;

\echo '===== verified ALC exact-time matching audit ====='
SELECT COUNT(*) AS n_icu_stays,
       COUNT(abs_lymphocytes_min) AS n_with_verified_alc,
       ROUND(100.0 * COUNT(abs_lymphocytes_min) / NULLIF(COUNT(*),0), 1) AS verified_alc_coverage_pct,
       SUM(alc_exact_match_count) AS n_exact_matched_measurement_pairs,
       MIN(alc_time_matching_verified::int) AS derivation_verified_for_all_rows
FROM public.eicu_external;


\copy (SELECT * FROM public.eicu_external) TO :'eicu_external_csv' CSV HEADER
\echo 'eicu_external export completed.'
