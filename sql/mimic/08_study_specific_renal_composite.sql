-- ============================================================
-- project_sa_aki : study-specific 30-day renal composite
-- This endpoint must not be described as a canonical MAKE30 definition.
-- Composite = death within 30 days OR new RRT OR persistent renal dysfunction.







-- ============================================================

SET search_path TO project_sa_aki, mimiciv_derived;

DROP TABLE IF EXISTS project_sa_aki.make30;

CREATE TABLE project_sa_aki.make30 AS
WITH

base AS (
    SELECT c.subject_id, c.hadm_id, c.stay_id, c.sepsis_onset_time,
           b.dischtime,
           COALESCE(o.mortality_30d, 0) AS mortality_30d
    FROM project_sa_aki.aki_kdigo_final c
    LEFT JOIN project_sa_aki.adult_icu_first b ON b.stay_id = c.stay_id
    LEFT JOIN project_sa_aki.outcomes o        ON o.stay_id = c.stay_id
),

rrt30 AS (
    SELECT c.stay_id,
           MAX(CASE WHEN r.charttime >= c.sepsis_onset_time
                     AND r.charttime <= c.sepsis_onset_time + INTERVAL '30 day'
                     AND (r.dialysis_present = 1 OR r.dialysis_active = 1)
                    THEN 1 ELSE 0 END) AS rrt_30d
    FROM project_sa_aki.aki_kdigo_final c
    LEFT JOIN mimiciv_derived.rrt r ON r.stay_id = c.stay_id
    GROUP BY c.stay_id
),

scr_window AS (
    SELECT stay_id, terminal_scr FROM (
        SELECT c.stay_id, ch.creatinine AS terminal_scr,
               ROW_NUMBER() OVER (PARTITION BY c.stay_id
                   ORDER BY ABS(EXTRACT(EPOCH FROM
                       (ch.charttime - (c.sepsis_onset_time + INTERVAL '30 day'))))) AS rn
        FROM base c
        JOIN mimiciv_derived.chemistry ch
          ON ch.subject_id = c.subject_id AND ch.creatinine IS NOT NULL
         AND ch.charttime >= c.sepsis_onset_time + INTERVAL '23 day'
         AND ch.charttime <= c.sepsis_onset_time + INTERVAL '37 day'
    ) t WHERE rn = 1
),

scr_disch AS (
    SELECT stay_id, terminal_scr FROM (
        SELECT c.stay_id, ch.creatinine AS terminal_scr,
               ROW_NUMBER() OVER (PARTITION BY c.stay_id
                   ORDER BY ch.charttime DESC) AS rn
        FROM base c
        JOIN mimiciv_derived.chemistry ch
          ON ch.subject_id = c.subject_id AND ch.creatinine IS NOT NULL
         AND ch.charttime >= c.sepsis_onset_time
         AND ch.charttime <= c.dischtime
        WHERE c.mortality_30d = 0
    ) t WHERE rn = 1
),

term AS (
    SELECT b.stay_id,
           COALESCE(sw.terminal_scr, sd.terminal_scr) AS terminal_scr,
           CASE WHEN sw.terminal_scr IS NOT NULL THEN 'window'
                WHEN sd.terminal_scr IS NOT NULL THEN 'discharge'
                ELSE NULL END AS terminal_scr_source
    FROM base b
    LEFT JOIN scr_window sw ON sw.stay_id = b.stay_id
    LEFT JOIN scr_disch  sd ON sd.stay_id = b.stay_id
),
prd AS (
    SELECT b.stay_id, t.terminal_scr, t.terminal_scr_source, cb.scr_baseline,
           CASE WHEN t.terminal_scr IS NOT NULL THEN 1 ELSE 0 END AS terminal_scr_available,
           CASE WHEN t.terminal_scr IS NOT NULL AND cb.scr_baseline IS NOT NULL
                 AND t.terminal_scr >= 1.5 * cb.scr_baseline
                THEN 1 ELSE 0 END AS persistent_rd_30d
    FROM base b
    LEFT JOIN term t ON t.stay_id = b.stay_id
    LEFT JOIN mimiciv_derived.creatinine_baseline cb ON cb.hadm_id = b.hadm_id
)
SELECT
    b.subject_id, b.hadm_id, b.stay_id,
    b.mortality_30d,
    COALESCE(r.rrt_30d, 0) AS rrt_30d,
    p.persistent_rd_30d,
    p.terminal_scr, p.scr_baseline,
    p.terminal_scr_available, p.terminal_scr_source,
    GREATEST(b.mortality_30d, COALESCE(r.rrt_30d,0), COALESCE(p.persistent_rd_30d,0)) AS make30
FROM base b
LEFT JOIN rrt30 r ON r.stay_id = b.stay_id
LEFT JOIN prd   p ON p.stay_id = b.stay_id;

CREATE INDEX IF NOT EXISTS idx_make30_stay ON project_sa_aki.make30 (stay_id);


SELECT COUNT(*) AS n, SUM(mortality_30d) AS n_death, SUM(rrt_30d) AS n_rrt,
       SUM(persistent_rd_30d) AS n_prd, SUM(make30) AS n_make30,
       SUM(terminal_scr_available) AS n_term_scr,
       ROUND(100.0*AVG(make30),1) AS make30_pct
FROM project_sa_aki.make30;

SELECT COALESCE(terminal_scr_source,'none') AS src, COUNT(*) AS n
FROM project_sa_aki.make30 GROUP BY 1 ORDER BY 2 DESC;
