-- ============================================================

\if :{?final_sepsis_csv}
\else
\set final_sepsis_csv './outputs/sql/final_sepsis_analysis.csv'
\endif




-- ============================================================

SET search_path TO project_sa_aki;

DROP TABLE IF EXISTS project_sa_aki.final_sepsis_analysis;

CREATE TABLE project_sa_aki.final_sepsis_analysis AS
SELECT

    base.subject_id, base.hadm_id, base.stay_id,


    base.gender, base.age, base.race, base.insurance,
    base.marital_status, base.admission_type,
    base.first_careunit, base.intime AS icu_intime,


    s3.sepsis_onset_time,
    s3.sofa_score,
    s3.respiration, s3.coagulation, s3.liver,
    s3.cardiovascular, s3.cns, s3.renal,


    aki.aki_creatinine_7d, aki.aki_stage_creatinine_7d,
    aki.aki_kdigo_7d,      aki.aki_stage_kdigo_7d,
    aki.rrt_24h, aki.rrt_48h, aki.rrt_7d,


    o.icu_mortality, o.hospital_mortality,
    o.mortality_28d, o.mortality_30d, o.mortality_90d,
    o.days_to_death,
    o.icu_los_days, o.hospital_los_days

FROM project_sa_aki.aki_kdigo_final aki
JOIN project_sa_aki.strict_sepsis3_cohort s3 ON s3.stay_id = aki.stay_id
JOIN project_sa_aki.adult_icu_first       base ON base.stay_id = aki.stay_id
LEFT JOIN project_sa_aki.outcomes         o    ON o.stay_id   = aki.stay_id;

CREATE INDEX IF NOT EXISTS idx_final_analysis_stay
    ON project_sa_aki.final_sepsis_analysis (stay_id);


SELECT
    COUNT(*) AS n,
    ROUND(100.0*AVG(aki_kdigo_7d),1)      AS aki_pct,
    ROUND(100.0*AVG(hospital_mortality),1) AS hosp_death_pct,
    ROUND(100.0*AVG(mortality_28d),1)     AS death28_pct,
    ROUND(AVG(sofa_score),1)              AS mean_sofa
FROM project_sa_aki.final_sepsis_analysis;



SELECT
    (SELECT COUNT(*) FROM project_sa_aki.aki_kdigo_final)        AS n_main_table,
    (SELECT COUNT(*) FROM project_sa_aki.final_sepsis_analysis)  AS n_wide_table,
    CASE WHEN (SELECT COUNT(*) FROM project_sa_aki.aki_kdigo_final)
            = (SELECT COUNT(*) FROM project_sa_aki.final_sepsis_analysis)
         THEN 'OK: matching row counts'
         ELSE 'WARNING: row count mismatch; check upstream tables for duplicated stay_id' END AS check_result;

-- ============================================================




-- ============================================================
\copy (SELECT * FROM project_sa_aki.final_sepsis_analysis) TO :'final_sepsis_csv' WITH CSV HEADER
