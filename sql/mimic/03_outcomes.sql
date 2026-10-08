-- ============================================================




-- LOS: icu_los_days=los; hospital_los_days=dischtime-admittime

-- ============================================================

SET search_path TO project_sa_aki, mimiciv_hosp, mimiciv_icu;

DROP TABLE IF EXISTS project_sa_aki.outcomes;

CREATE TABLE project_sa_aki.outcomes AS
SELECT
    c.subject_id, c.hadm_id, c.stay_id,


    CASE WHEN c.dod IS NOT NULL
          AND c.dod >= c.intime::date
          AND c.dod <= c.outtime::date
         THEN 1 ELSE 0 END                             AS icu_mortality,


    c.hospital_expire_flag::int                        AS hospital_mortality,


    CASE WHEN c.dod IS NOT NULL
          AND (c.dod - c.admittime::date) <= 28 THEN 1 ELSE 0 END AS mortality_28d,
    CASE WHEN c.dod IS NOT NULL
          AND (c.dod - c.admittime::date) <= 30 THEN 1 ELSE 0 END AS mortality_30d,
    CASE WHEN c.dod IS NOT NULL
          AND (c.dod - c.admittime::date) <= 90 THEN 1 ELSE 0 END AS mortality_90d,


    CASE WHEN c.dod IS NOT NULL
         THEN (c.dod - c.admittime::date) END          AS days_to_death,


    ROUND(c.los::numeric, 2)                           AS icu_los_days,
    ROUND(EXTRACT(EPOCH FROM (c.dischtime - c.admittime))/86400.0, 2)
                                                       AS hospital_los_days
FROM project_sa_aki.adult_icu_first c;

CREATE INDEX IF NOT EXISTS idx_outcomes_stay
    ON project_sa_aki.outcomes (stay_id);


SELECT
    COUNT(*) AS n,
    SUM(icu_mortality)     AS n_icu_death,
    SUM(hospital_mortality) AS n_hosp_death,
    SUM(mortality_28d)     AS n_28d,
    SUM(mortality_30d)     AS n_30d,
    SUM(mortality_90d)     AS n_90d,
    ROUND(AVG(icu_los_days),2)      AS mean_icu_los,
    ROUND(AVG(hospital_los_days),2) AS mean_hosp_los
FROM project_sa_aki.outcomes;
