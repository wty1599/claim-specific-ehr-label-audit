-- ============================================================
-- project_sa_aki : assemble outcomes + 24-hour features
-- ============================================================
-- Recovered production script: merge_export_full.sql
-- Provenance limitation: the original standalone producer was not present
-- in the audited source roots; this block was recovered from the prior
-- submission-code release and therefore requires author verification.
-- This creates the initial 80-column final_full table. The later
-- renal-composite merge adds six columns and produces the locked
-- 20,049 x 86 analysis CSV.
-- ============================================================

SET search_path TO project_sa_aki;

DROP TABLE IF EXISTS project_sa_aki.final_full;

CREATE TABLE project_sa_aki.final_full AS
SELECT
    f.*,
    feat.wbc_max, feat.wbc_min, feat.platelets_min, feat.hemoglobin_min,
    feat.abs_neutrophils_max, feat.abs_lymphocytes_min, feat.bands_max,
    feat.bilirubin_total_max, feat.alt_max, feat.ast_max, feat.albumin_min,
    feat.inr_max, feat.pt_max, feat.ptt_max, feat.fibrinogen_min,
    feat.creatinine_max, feat.creatinine_min, feat.bun_max, feat.aniongap_max,
    feat.bicarbonate_min, feat.glucose_max,
    feat.sodium_min, feat.sodium_max, feat.potassium_min, feat.potassium_max,
    feat.calcium_min, feat.calcium_max, feat.chloride_min, feat.chloride_max,
    feat.lactate_max, feat.ph_min, feat.po2_min, feat.pco2_max, feat.so2_min,
    feat.pao2fio2ratio_min, feat.baseexcess_min,
    feat.heart_rate_max, feat.sbp_min, feat.mbp_min, feat.mbp_mean,
    feat.resp_rate_max, feat.temperature_min, feat.temperature_max, feat.spo2_min,
    feat.gcs_min, feat.urine_output_24h_ml
FROM project_sa_aki.final_sepsis_analysis f
LEFT JOIN project_sa_aki.features_24h feat ON feat.stay_id = f.stay_id;

CREATE INDEX IF NOT EXISTS idx_final_full_stay
    ON project_sa_aki.final_full (stay_id);

SELECT
    (SELECT COUNT(*) FROM project_sa_aki.final_sepsis_analysis) AS n_base,
    (SELECT COUNT(*) FROM project_sa_aki.final_full) AS n_full,
    (SELECT COUNT(*) FROM information_schema.columns
      WHERE table_schema = 'project_sa_aki'
        AND table_name = 'final_full') AS n_cols;

-- Export is intentionally left to the local psql environment:
-- \copy (SELECT * FROM project_sa_aki.final_full)
--   TO '/path/to/output/final_full.csv' WITH CSV HEADER
