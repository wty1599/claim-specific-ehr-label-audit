-- ============================================================
-- build_baseline_covars_v2.sql




-- ============================================================
\set ON_ERROR_STOP on
\if :{?baseline_covars_csv}
\else
\set baseline_covars_csv './outputs/sql/baseline_covars.csv'
\endif

DROP TABLE IF EXISTS project_sa_aki.baseline_covars;

CREATE TABLE project_sa_aki.baseline_covars AS
WITH aki_24h AS (
  SELECT ie.stay_id,
         MAX(k.aki_stage_creat) AS stage_creat_0_24h,
         MAX(k.aki_stage_uo)    AS stage_uo_0_24h
  FROM mimiciv_icu.icustays ie
  JOIN mimiciv_derived.kdigo_stages k
    ON k.stay_id   = ie.stay_id
   AND k.charttime >= ie.intime
   AND k.charttime <  ie.intime + INTERVAL '24 hour'
  GROUP BY ie.stay_id
),
rrt_24h AS (
  SELECT ie.stay_id, 1 AS rrt_flag
  FROM mimiciv_icu.icustays ie
  JOIN mimiciv_derived.rrt r
    ON r.stay_id   = ie.stay_id
   AND r.charttime >= ie.intime
   AND r.charttime <  ie.intime + INTERVAL '24 hour'
   AND r.dialysis_active = 1
  GROUP BY ie.stay_id
)
SELECT ie.stay_id,
       p.gender,
       c.charlson_comorbidity_index,
       s.sapsii,
       GREATEST(
         COALESCE(a.stage_creat_0_24h, 0),
         COALESCE(a.stage_uo_0_24h,    0),
         CASE WHEN rr.rrt_flag = 1 THEN 3 ELSE 0 END
       ) AS aki_stage_0_24h
FROM mimiciv_icu.icustays ie
JOIN      mimiciv_hosp.patients   p ON p.subject_id = ie.subject_id
LEFT JOIN mimiciv_derived.charlson c ON c.hadm_id   = ie.hadm_id
LEFT JOIN mimiciv_derived.sapsii   s ON s.stay_id   = ie.stay_id
LEFT JOIN aki_24h  a  ON a.stay_id  = ie.stay_id
LEFT JOIN rrt_24h  rr ON rr.stay_id = ie.stay_id;

-- ---- sanity checks ----
\echo '
SELECT aki_stage_0_24h, COUNT(*) AS n
FROM project_sa_aki.baseline_covars GROUP BY 1 ORDER BY 1;

\echo '--- Completeness (nonmissing counts per column, expected near n) ---'
SELECT COUNT(*) AS n,
       COUNT(gender) AS n_gender,
       COUNT(charlson_comorbidity_index) AS n_charlson,
       COUNT(sapsii) AS n_sapsii
FROM project_sa_aki.baseline_covars;

\echo '
SELECT MIN(sapsii) AS min, ROUND(AVG(sapsii)) AS mean,
       PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY sapsii) AS median,
       MAX(sapsii) AS max
FROM project_sa_aki.baseline_covars WHERE sapsii IS NOT NULL;


\copy (SELECT * FROM project_sa_aki.baseline_covars) TO :'baseline_covars_csv' CSV HEADER
\echo 'baseline_covars export completed.'
