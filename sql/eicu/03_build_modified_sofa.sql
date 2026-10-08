-- ============================================================










-- ============================================================
\set ON_ERROR_STOP on
\pset pager off
\if :{?eicu_external_csv}
\else
\set eicu_external_csv './outputs/sql/eicu_external.csv'
\endif


DROP TABLE IF EXISTS tmp_plt;
CREATE TEMP TABLE tmp_plt AS
SELECT l.patientunitstayid, MIN(l.labresult) AS plt_min
FROM lab l JOIN public.saaki_cohort c ON c.patientunitstayid=l.patientunitstayid
WHERE lower(l.labname)='platelets x 1000' AND l.labresultoffset BETWEEN 0 AND 1440
  AND l.labresult>0
GROUP BY l.patientunitstayid;


DROP TABLE IF EXISTS tmp_sf;
CREATE TEMP TABLE tmp_sf AS
SELECT c.patientunitstayid,
  (SELECT MIN(v.sao2) FROM vitalperiodic v
    WHERE v.patientunitstayid=c.patientunitstayid AND v.observationoffset BETWEEN 0 AND 1440
      AND v.sao2 BETWEEN 50 AND 100) AS sao2_min,
  (SELECT MAX(l.labresult) FROM lab l
    WHERE l.patientunitstayid=c.patientunitstayid AND lower(l.labname)='fio2'
      AND l.labresultoffset BETWEEN 0 AND 1440 AND l.labresult BETWEEN 21 AND 100) AS fio2_max
FROM public.saaki_cohort c;


DROP TABLE IF EXISTS public.eicu_sofa;
CREATE TABLE public.eicu_sofa AS
WITH ap AS (
  SELECT patientunitstayid,
    NULLIF(meanbp,-1)     AS meanbp,
    NULLIF(bilirubin,-1)  AS bili,
    NULLIF(creatinine,-1) AS creat,
    CASE WHEN eyes IS NOT NULL AND motor IS NOT NULL AND verbal IS NOT NULL
         THEN eyes+motor+verbal END AS gcs
  FROM apacheapsvar
),
comp AS (
  SELECT c.patientunitstayid,

    CASE WHEN sf.sao2_min IS NOT NULL AND sf.fio2_max IS NOT NULL AND sf.fio2_max>0 THEN
      (CASE WHEN (sf.sao2_min/(sf.fio2_max/100.0)) <= 67  THEN 4
            WHEN (sf.sao2_min/(sf.fio2_max/100.0)) <= 142 THEN 3
            WHEN (sf.sao2_min/(sf.fio2_max/100.0)) <= 221 THEN 2
            WHEN (sf.sao2_min/(sf.fio2_max/100.0)) <= 301 THEN 1 ELSE 0 END)
    ELSE 0 END AS resp,

    CASE WHEN pl.plt_min IS NULL THEN 0
         WHEN pl.plt_min < 20  THEN 4 WHEN pl.plt_min < 50 THEN 3
         WHEN pl.plt_min < 100 THEN 2 WHEN pl.plt_min < 150 THEN 1 ELSE 0 END AS coag,

    CASE WHEN ap.bili IS NULL THEN 0
         WHEN ap.bili >= 12 THEN 4 WHEN ap.bili >= 6 THEN 3
         WHEN ap.bili >= 2  THEN 2 WHEN ap.bili >= 1.2 THEN 1 ELSE 0 END AS liver,

    CASE WHEN ap.meanbp IS NULL THEN 0 WHEN ap.meanbp < 70 THEN 1 ELSE 0 END AS cardio,

    CASE WHEN ap.gcs IS NULL THEN 0
         WHEN ap.gcs < 6 THEN 4 WHEN ap.gcs < 10 THEN 3
         WHEN ap.gcs < 13 THEN 2 WHEN ap.gcs < 15 THEN 1 ELSE 0 END AS cns,

    CASE WHEN ap.creat IS NULL THEN 0
         WHEN ap.creat >= 5 THEN 4 WHEN ap.creat >= 3.5 THEN 3
         WHEN ap.creat >= 2 THEN 2 WHEN ap.creat >= 1.2 THEN 1 ELSE 0 END AS renal
  FROM public.saaki_cohort c
  LEFT JOIN ap    ON ap.patientunitstayid=c.patientunitstayid
  LEFT JOIN tmp_plt pl ON pl.patientunitstayid=c.patientunitstayid
  LEFT JOIN tmp_sf sf  ON sf.patientunitstayid=c.patientunitstayid
)
SELECT patientunitstayid, resp, coag, liver, cardio, cns, renal,
       (resp+coag+liver+cardio+cns+renal) AS sofa_total
FROM comp;

\echo '===== Modified SOFA distribution ====='
SELECT ROUND(MIN(sofa_total),0) min, ROUND(AVG(sofa_total),1) mean,
       PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY sofa_total) median,
       ROUND(MAX(sofa_total),0) max FROM public.eicu_sofa;
\echo '===== Component means (plausibility check) ====='
SELECT ROUND(AVG(resp),2) resp, ROUND(AVG(coag),2) coag, ROUND(AVG(liver),2) liver,
       ROUND(AVG(cardio),2) cardio, ROUND(AVG(cns),2) cns, ROUND(AVG(renal),2) renal
FROM public.eicu_sofa;


ALTER TABLE public.eicu_external DROP COLUMN IF EXISTS sofa;
ALTER TABLE public.eicu_external ADD COLUMN sofa numeric;
UPDATE public.eicu_external e SET sofa = s.sofa_total
FROM public.eicu_sofa s WHERE s.patientunitstayid = e.patientunitstayid;

\echo '===== Re-export eicu_external(including sofa) ====='
\copy (SELECT * FROM public.eicu_external) TO :'eicu_external_csv' CSV HEADER
\echo 'Export updated'
