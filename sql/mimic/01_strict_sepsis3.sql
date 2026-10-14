-- ============================================================







-- ============================================================

SET search_path TO project_sa_aki, mimiciv_derived, mimiciv_icu, mimiciv_hosp;

DROP TABLE IF EXISTS project_sa_aki.strict_sepsis3_cohort;

CREATE TABLE project_sa_aki.strict_sepsis3_cohort AS
WITH sofa_ge2 AS (
    SELECT stay_id, starttime, endtime,
           respiration_24hours  AS respiration,
           coagulation_24hours  AS coagulation,
           liver_24hours        AS liver,
           cardiovascular_24hours AS cardiovascular,
           cns_24hours          AS cns,
           renal_24hours        AS renal,
           sofa_24hours         AS sofa_score
    FROM mimiciv_derived.sofa
    WHERE sofa_24hours >= 2
),
s1 AS (
    SELECT
        soi.subject_id, soi.stay_id,
        soi.antibiotic_time, soi.culture_time,
        soi.suspected_infection_time,
        sf.endtime AS sofa_time,
        sf.sofa_score,
        sf.respiration, sf.coagulation, sf.liver,
        sf.cardiovascular, sf.cns, sf.renal,
        ROW_NUMBER() OVER (
            PARTITION BY soi.stay_id
            ORDER BY soi.suspected_infection_time, soi.antibiotic_time,
                     soi.culture_time, sf.endtime
        ) AS rn_sus
    FROM mimiciv_derived.suspicion_of_infection soi
    JOIN sofa_ge2 sf
      ON sf.stay_id = soi.stay_id
     AND sf.endtime >= soi.suspected_infection_time - INTERVAL '48 hour'   -- (b)
     AND sf.endtime <= soi.suspected_infection_time + INTERVAL '24 hour'
    WHERE soi.stay_id IS NOT NULL
      AND soi.suspected_infection = 1
)
SELECT
    c.*,
    s1.suspected_infection_time,
    s1.suspected_infection_time AS sepsis_onset_time,
    s1.antibiotic_time, s1.culture_time,
    s1.sofa_time, s1.sofa_score,
    s1.respiration, s1.coagulation, s1.liver,
    s1.cardiovascular, s1.cns, s1.renal,
    1::int AS sepsis3
FROM project_sa_aki.adult_icu_first c
JOIN s1 ON s1.stay_id = c.stay_id
WHERE s1.rn_sus = 1
  AND s1.suspected_infection_time
        BETWEEN c.icu_intime - INTERVAL '24 hour'          -- Cohort window: 24 hours before ICU admission.
            AND c.icu_intime + INTERVAL '24 hour';

CREATE INDEX IF NOT EXISTS idx_strict_sepsis3_stay
    ON project_sa_aki.strict_sepsis3_cohort (stay_id);


SELECT COUNT(*) AS n_sepsis3,
       MIN(sofa_score) AS sofa_min,
       ROUND(AVG(sofa_score),2) AS sofa_mean,
       MAX(sofa_score) AS sofa_max
FROM project_sa_aki.strict_sepsis3_cohort;


-- ============================================================
-- Optional diagnostic query to quantify exclusions from the cohort definition.
-- ============================================================




-- WITH cohort_soi AS (
--     SELECT c.stay_id, c.subject_id, c.icu_intime,
--            soi.suspected_infection_time AS si_time,
--            soi.stay_id AS soi_stay_id
--     FROM project_sa_aki.adult_icu_first c
--     JOIN mimiciv_derived.suspicion_of_infection soi
--       ON soi.subject_id = c.subject_id
--      AND soi.suspected_infection = 1
-- )
-- SELECT
--     COUNT(*) FILTER (WHERE soi_stay_id IS NULL
--                       AND si_time < icu_intime)      AS n_pre_icu_soi_excluded,
--     COUNT(*) FILTER (WHERE soi_stay_id IS NOT NULL)  AS n_in_icu_soi,
--     COUNT(*)                                          AS n_total_soi_rows
-- FROM cohort_soi;


-- SELECT
--   (SELECT COUNT(*) FROM project_sa_aki.adult_icu_first)            AS step0_adult_icu_first,
--   (SELECT COUNT(DISTINCT stay_id) FROM mimiciv_derived.suspicion_of_infection
--      WHERE suspected_infection = 1 AND stay_id IS NOT NULL)         AS step1_with_soi_in_icu,
--   (SELECT COUNT(*) FROM project_sa_aki.strict_sepsis3_cohort)      AS step2_strict_sepsis3;


-- SELECT
--   ROUND(AVG(respiration),2)   AS resp,
--   ROUND(AVG(coagulation),2)   AS coag,
--   ROUND(AVG(liver),2)         AS liver,
--   ROUND(AVG(cardiovascular),2) AS cardio,
--   ROUND(AVG(cns),2)           AS cns,
--   ROUND(AVG(renal),2)         AS renal
-- FROM project_sa_aki.strict_sepsis3_cohort;
