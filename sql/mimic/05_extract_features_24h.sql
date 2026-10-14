-- ============================================================





-- ============================================================

SET search_path TO project_sa_aki, mimiciv_derived;

DROP TABLE IF EXISTS project_sa_aki.features_24h;

CREATE TABLE project_sa_aki.features_24h AS
SELECT
    c.subject_id, c.hadm_id, c.stay_id,


    lab.wbc_max, lab.wbc_min,
    lab.platelets_min,
    lab.hemoglobin_min,
    lab.abs_neutrophils_max, lab.abs_lymphocytes_min, lab.bands_max,


    lab.bilirubin_total_max,
    lab.alt_max, lab.ast_max,
    lab.albumin_min,


    lab.inr_max, lab.pt_max, lab.ptt_max, lab.fibrinogen_min,


    lab.creatinine_max, lab.creatinine_min,
    lab.bun_max,
    lab.aniongap_max,
    lab.bicarbonate_min,
    lab.glucose_max,

    lab.sodium_min, lab.sodium_max,
    lab.potassium_min, lab.potassium_max,
    lab.calcium_min, lab.calcium_max,
    lab.chloride_min, lab.chloride_max,


    bg.lactate_max,
    bg.ph_min,
    bg.po2_min, bg.pco2_max, bg.so2_min,
    bg.pao2fio2ratio_min,
    bg.baseexcess_min,


    vs.heart_rate_max,
    vs.sbp_min, vs.mbp_min, vs.mbp_mean,
    vs.resp_rate_max,
    vs.temperature_min, vs.temperature_max,
    vs.spo2_min,


    gcs.gcs_min,


    uo.urineoutput AS urine_output_24h_ml

FROM project_sa_aki.final_sepsis_analysis c
LEFT JOIN mimiciv_derived.first_day_lab          lab ON lab.stay_id = c.stay_id
LEFT JOIN mimiciv_derived.first_day_bg           bg  ON bg.stay_id  = c.stay_id
LEFT JOIN mimiciv_derived.first_day_vitalsign    vs  ON vs.stay_id  = c.stay_id
LEFT JOIN mimiciv_derived.first_day_gcs          gcs ON gcs.stay_id = c.stay_id
LEFT JOIN mimiciv_derived.first_day_urine_output uo  ON uo.stay_id  = c.stay_id;

CREATE INDEX IF NOT EXISTS idx_features24h_stay
    ON project_sa_aki.features_24h (stay_id);


SELECT
    (SELECT COUNT(*) FROM project_sa_aki.final_sepsis_analysis) AS n_main,
    (SELECT COUNT(*) FROM project_sa_aki.features_24h)          AS n_feat,
    CASE WHEN (SELECT COUNT(*) FROM project_sa_aki.final_sepsis_analysis)
            = (SELECT COUNT(*) FROM project_sa_aki.features_24h)
         THEN 'OK' ELSE 'WARN row count mismatch' END AS check;

-- Missingness percentages for selected variables.
SELECT
    ROUND(100.0*SUM(CASE WHEN lactate_max      IS NULL THEN 1 ELSE 0 END)/COUNT(*),1) AS miss_lactate,
    ROUND(100.0*SUM(CASE WHEN creatinine_max   IS NULL THEN 1 ELSE 0 END)/COUNT(*),1) AS miss_creat,
    ROUND(100.0*SUM(CASE WHEN gcs_min          IS NULL THEN 1 ELSE 0 END)/COUNT(*),1) AS miss_gcs,
    ROUND(100.0*SUM(CASE WHEN bands_max        IS NULL THEN 1 ELSE 0 END)/COUNT(*),1) AS miss_bands,
    ROUND(100.0*SUM(CASE WHEN albumin_min      IS NULL THEN 1 ELSE 0 END)/COUNT(*),1) AS miss_albumin,
    ROUND(100.0*SUM(CASE WHEN pao2fio2ratio_min IS NULL THEN 1 ELSE 0 END)/COUNT(*),1) AS miss_pf,
    ROUND(100.0*SUM(CASE WHEN urine_output_24h_ml IS NULL THEN 1 ELSE 0 END)/COUNT(*),1) AS miss_uo
FROM project_sa_aki.features_24h;
