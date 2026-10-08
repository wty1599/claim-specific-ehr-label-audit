-- ============================================================

-- ============================================================


-- Public-release terminology: eICU sepsis is a diagnosis-string proxy rather
-- than a direct reimplementation of the MIMIC Sepsis-3 onset algorithm.
\set ON_ERROR_STOP on
\pset pager off

DROP TABLE IF EXISTS public.saaki_cohort;

CREATE TEMP TABLE tmp_sepsis AS
SELECT DISTINCT patientunitstayid FROM diagnosis
WHERE lower(diagnosisstring) LIKE '%sepsis%' OR lower(diagnosisstring) LIKE '%septic%';



CREATE TEMP TABLE tmp_ckd_esrd AS
SELECT DISTINCT patientunitstayid FROM pasthistory
WHERE lower(pasthistoryvalue) IN (
    'renal failure - hemodialysis',
    'renal failure - peritoneal dialysis',
    'renal failure- not currently dialyzed',
    's/p renal transplant',
    'renal insufficiency - creatinine 2-3',
    'renal insufficiency - creatinine 3-4',
    'renal insufficiency - creatinine 4-5',
    'renal insufficiency - creatinine > 5'
)
UNION
SELECT DISTINCT patientunitstayid FROM diagnosis
WHERE lower(diagnosisstring) LIKE 'renal|disorder of kidney|esrd%'
   OR lower(diagnosisstring) LIKE 'renal|disorder of kidney|chronic kidney disease%'
   OR lower(diagnosisstring) LIKE 'renal|disorder of kidney|chronic renal insufficiency%'
   OR icd9code LIKE '%585.6%' OR icd9code LIKE '%V45.11%'
   OR icd9code LIKE '%N18.5%' OR icd9code LIKE '%N18.6%' OR icd9code LIKE '%Z99.2%';

CREATE TEMP TABLE tmp_scr_meas AS
SELECT patientunitstayid, MIN(labresult) AS scr_meas_low
FROM lab
WHERE lower(labname)='creatinine' AND labresultoffset BETWEEN -360 AND 1440
  AND labresult>0 AND labresult<30
GROUP BY patientunitstayid;

CREATE TEMP TABLE tmp_scr_peak AS
SELECT patientunitstayid, MAX(labresult) AS scr_peak_0_24h
FROM lab
WHERE lower(labname)='creatinine' AND labresultoffset BETWEEN 0 AND 1440
  AND labresult>0 AND labresult<30
GROUP BY patientunitstayid;

CREATE TEMP TABLE tmp_rrt_24h AS
SELECT DISTINCT patientunitstayid FROM treatment
WHERE (lower(treatmentstring) LIKE '%dialysis%' OR lower(treatmentstring) LIKE '%crrt%'
    OR lower(treatmentstring) LIKE '%c v v h%' OR lower(treatmentstring) LIKE '%hemodialysis%'
    OR lower(treatmentstring) LIKE '%ultrafiltration%' OR lower(treatmentstring) LIKE '%renal replacement%')
  AND treatmentoffset BETWEEN 0 AND 1440;


CREATE TEMP TABLE tmp_demo AS
SELECT p.patientunitstayid,
  CASE WHEN p.age='> 89' THEN 90 WHEN p.age ~ '^[0-9]+$' THEN p.age::numeric END AS age_n,
  CASE WHEN lower(p.gender)='female' THEN 1 ELSE 0 END AS is_female,
  CASE WHEN lower(p.ethnicity) LIKE '%african%' OR lower(p.ethnicity) LIKE '%black%' THEN 1 ELSE 0 END AS is_black
FROM patient p;

CREATE TEMP TABLE tmp_mdrd AS
SELECT patientunitstayid,
  POWER( 75.0 / (POWER(1.212,is_black) * POWER(0.742,is_female)) / 186.0
         * POWER(NULLIF(age_n,0), 0.203), -1.0/1.154 ) AS scr_mdrd
FROM tmp_demo WHERE age_n IS NOT NULL;


CREATE TEMP TABLE tmp_aki2 AS
SELECT pk.patientunitstayid, pk.scr_peak_0_24h,
       LEAST(COALESCE(m.scr_mdrd, ms.scr_meas_low), COALESCE(ms.scr_meas_low, m.scr_mdrd)) AS scr_baseline,
       CASE WHEN r.patientunitstayid IS NOT NULL THEN 1 ELSE 0 END AS rrt_flag
FROM tmp_scr_peak pk
LEFT JOIN tmp_scr_meas ms ON ms.patientunitstayid=pk.patientunitstayid
LEFT JOIN tmp_mdrd m      ON m.patientunitstayid=pk.patientunitstayid
LEFT JOIN tmp_rrt_24h r   ON r.patientunitstayid=pk.patientunitstayid;

CREATE TABLE public.saaki_cohort AS
SELECT
  p.patientunitstayid, p.uniquepid, p.hospitalid,
  CASE WHEN p.age='> 89' THEN 90 WHEN p.age ~ '^[0-9]+$' THEN p.age::int END AS age,
  CASE WHEN lower(p.gender)='male' THEN 1 WHEN lower(p.gender)='female' THEN 0 END AS sex,
  ap.apachescore AS apache_iva,
  a.scr_baseline, a.scr_peak_0_24h,
  CASE
    WHEN a.rrt_flag=1 THEN 3
    WHEN a.scr_peak_0_24h >= 3.0*a.scr_baseline OR a.scr_peak_0_24h >= 4.0 THEN 3
    WHEN a.scr_peak_0_24h >= 2.0*a.scr_baseline THEN 2
    WHEN a.scr_peak_0_24h >= 1.5*a.scr_baseline OR a.scr_peak_0_24h >= a.scr_baseline+0.3 THEN 1
    ELSE 0 END AS aki_stage_0_24h,
  CASE WHEN lower(p.hospitaldischargestatus)='expired' THEN 1
       WHEN lower(p.hospitaldischargestatus)='alive' THEN 0 END AS hosp_mortality,
  CASE WHEN lower(p.unitdischargestatus)='expired' THEN 1 ELSE 0 END AS icu_mortality,
  p.hospitaldischargeoffset, p.unitdischargeoffset,
  dc.scr_discharge
FROM patient p
JOIN tmp_sepsis s ON s.patientunitstayid=p.patientunitstayid
JOIN tmp_aki2 a   ON a.patientunitstayid=p.patientunitstayid
LEFT JOIN apachepatientresult ap ON ap.patientunitstayid=p.patientunitstayid AND ap.apachescore>=0
LEFT JOIN LATERAL (
  SELECT labresult AS scr_discharge FROM lab l
  WHERE l.patientunitstayid=p.patientunitstayid AND lower(l.labname)='creatinine'
    AND l.labresult>0 AND l.labresult<30
  ORDER BY l.labresultoffset DESC LIMIT 1) dc ON true
WHERE (p.age='> 89' OR (p.age ~ '^[0-9]+$' AND p.age::int>=18))
  AND p.unitvisitnumber=1
  AND p.patientunitstayid NOT IN (SELECT patientunitstayid FROM tmp_ckd_esrd);

\echo '===== [1] Exclusion flow counts ====='
SELECT
  (SELECT COUNT(DISTINCT patientunitstayid) FROM tmp_sepsis)                          AS a_sepsis,
  (SELECT COUNT(*) FROM tmp_ckd_esrd)                                                  AS b_ckd_esrd_flagged,
  (SELECT COUNT(*) FROM public.saaki_cohort)                                           AS c_final_cohort,
  (SELECT COUNT(DISTINCT hospitalid) FROM public.saaki_cohort)                         AS n_hospitals;

\echo '===== [2] AKI stage distribution (historical expectation: stage 0 populated, not 0%) ====='
SELECT aki_stage_0_24h, COUNT(*),
       ROUND(100.0*COUNT(*)/SUM(COUNT(*)) OVER (),1) AS pct
FROM public.saaki_cohort GROUP BY 1 ORDER BY 1;

\echo '===== [3] In-hospital mortality (historical expectation: below ~0.236 after retaining sepsis without AKI) ====='
SELECT ROUND(AVG(hosp_mortality)::numeric,3) AS hosp_mort,
       COUNT(*) AS n FROM public.saaki_cohort;

\echo '===== [4] Baseline creatinine distribution (historical expectation: median ~0.8-1.2; MDRD back-calculation unchanged) ====='
SELECT ROUND((PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY scr_baseline))::numeric,2) AS median_scr
FROM public.saaki_cohort;
SELECT COUNT(*) n, COUNT(apache_iva) n_apache, ROUND(AVG(hosp_mortality)::numeric,3) hosp_mort FROM public.saaki_cohort;
