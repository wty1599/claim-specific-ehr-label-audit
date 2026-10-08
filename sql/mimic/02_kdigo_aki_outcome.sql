-- ============================================================




--       rrt_24h / rrt_48h / rrt_7d


-- ============================================================

SET search_path TO project_sa_aki, mimiciv_derived, mimiciv_hosp, mimiciv_icu;

DROP TABLE IF EXISTS project_sa_aki.aki_kdigo_final;

CREATE TABLE project_sa_aki.aki_kdigo_final AS
WITH



esrd_dx AS (
    SELECT DISTINCT hadm_id
    FROM mimiciv_hosp.diagnoses_icd
    WHERE (icd_version = 9  AND icd_code IN ('5856','V4511'))
       OR (icd_version = 10 AND icd_code IN ('N186','N185','Z992'))
),
esrd_flag AS (
    SELECT c.stay_id, c.hadm_id,
           CASE WHEN cb.ckd = 1 OR e.hadm_id IS NOT NULL THEN 1 ELSE 0 END AS esrd_exclude
    FROM project_sa_aki.strict_sepsis3_cohort c
    LEFT JOIN mimiciv_derived.creatinine_baseline cb ON cb.hadm_id = c.hadm_id
    LEFT JOIN esrd_dx e ON e.hadm_id = c.hadm_id
),

ks_win AS (
    SELECT
        c.stay_id,
        ks.aki_stage_creat,
        ks.aki_stage_uo,
        ks.aki_stage_crrt
    FROM project_sa_aki.strict_sepsis3_cohort c
    JOIN mimiciv_derived.kdigo_stages ks
      ON ks.stay_id = c.stay_id
     AND ks.charttime >= c.sepsis_onset_time
     AND ks.charttime <= c.sepsis_onset_time + INTERVAL '7 day'
),

aki_agg AS (
    SELECT
        stay_id,
        MAX(aki_stage_creat) AS aki_stage_creatinine_7d,
        MAX(aki_stage_uo)    AS aki_stage_uo_7d,

        MAX(GREATEST(
              COALESCE(aki_stage_creat,0),
              COALESCE(aki_stage_uo,0),
              COALESCE(aki_stage_crrt,0)
        )) AS aki_stage_kdigo_7d
    FROM ks_win
    GROUP BY stay_id
),

rrt_win AS (
    SELECT
        c.stay_id,
        MAX(CASE WHEN r.charttime >= c.sepsis_onset_time
                  AND r.charttime <= c.sepsis_onset_time + INTERVAL '24 hour'
                  AND (r.dialysis_present = 1 OR r.dialysis_active = 1)
                 THEN 1 ELSE 0 END) AS rrt_24h,
        MAX(CASE WHEN r.charttime >= c.sepsis_onset_time
                  AND r.charttime <= c.sepsis_onset_time + INTERVAL '48 hour'
                  AND (r.dialysis_present = 1 OR r.dialysis_active = 1)
                 THEN 1 ELSE 0 END) AS rrt_48h,
        MAX(CASE WHEN r.charttime >= c.sepsis_onset_time
                  AND r.charttime <= c.sepsis_onset_time + INTERVAL '7 day'
                  AND (r.dialysis_present = 1 OR r.dialysis_active = 1)
                 THEN 1 ELSE 0 END) AS rrt_7d
    FROM project_sa_aki.strict_sepsis3_cohort c
    LEFT JOIN mimiciv_derived.rrt r ON r.stay_id = c.stay_id
    GROUP BY c.stay_id
)

SELECT
    c.subject_id, c.hadm_id, c.stay_id,
    c.sepsis_onset_time,

    COALESCE(a.aki_stage_creatinine_7d, 0)                    AS aki_stage_creatinine_7d,
    (COALESCE(a.aki_stage_creatinine_7d, 0) >= 1)::int        AS aki_creatinine_7d,

    COALESCE(a.aki_stage_kdigo_7d, 0)                         AS aki_stage_kdigo_7d,
    (COALESCE(a.aki_stage_kdigo_7d, 0) >= 1)::int             AS aki_kdigo_7d,

    COALESCE(rw.rrt_24h, 0) AS rrt_24h,
    COALESCE(rw.rrt_48h, 0) AS rrt_48h,
    COALESCE(rw.rrt_7d, 0)  AS rrt_7d
FROM project_sa_aki.strict_sepsis3_cohort c
JOIN esrd_flag ef ON ef.stay_id = c.stay_id
LEFT JOIN aki_agg  a  ON a.stay_id  = c.stay_id
LEFT JOIN rrt_win  rw ON rw.stay_id = c.stay_id
WHERE ef.esrd_exclude = 0;

CREATE INDEX IF NOT EXISTS idx_aki_kdigo_stay
    ON project_sa_aki.aki_kdigo_final (stay_id);


SELECT
    COUNT(*)                              AS n_cohort_after_esrd_excl,
    SUM(aki_creatinine_7d)                AS n_aki_creat,
    SUM(aki_kdigo_7d)                     AS n_aki_kdigo,
    SUM(rrt_7d)                           AS n_rrt_7d,
    ROUND(100.0*SUM(aki_kdigo_7d)/COUNT(*),1) AS aki_kdigo_pct
FROM project_sa_aki.aki_kdigo_final;

SELECT aki_stage_kdigo_7d, COUNT(*) AS n
FROM project_sa_aki.aki_kdigo_final
GROUP BY aki_stage_kdigo_7d ORDER BY aki_stage_kdigo_7d;
