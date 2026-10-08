\set ON_ERROR_STOP on

\if :{?rrt_first_csv}
\else
\set rrt_first_csv './outputs/sql/rrt_first.csv'
\endif
\if :{?rrt_windows_csv}
\else
\set rrt_windows_csv './outputs/sql/rrt_timing_windows_v2.csv'
\endif
\if :{?rrt_base_times_csv}
\else
\set rrt_base_times_csv './outputs/sql/rrt_base_times.csv'
\endif
\if :{?features_0_6_csv}
\else
\set features_0_6_csv './outputs/sql/features_0_6.csv'
\endif

-- =============================================================
-- Public release: output paths are psql variables. Create their parent
-- directories before execution or override all four variables with -v.
-- =============================================================


-- =============================================================
-- 00_extract_rrt_hte_features_abs.sql
-- Purpose:
--   Extract data required for revised RRT HTE analyses in EHR audit:
--   1) first RRT time relative to ICU admission
--   2) RRT time windows: before ICU, 0-6h, 6-24h, 24-72h, >72h, no RRT
--   3) death/landmark time variables
--   4) 0-6h features for early phenotype and 6-24h RRT HTE
--
-- Assumptions:
--   - Existing cohort table: project_sa_aki.final_full
--   - ICU stay table: mimiciv_icu.icustays
--   - Derived tables from mimic-code exist under mimiciv_derived:
--       rrt, chemistry, complete_blood_count, coagulation, enzyme,
--       bg, vitalsign, urine_output, gcs
--
-- Output CSV files:
--   ${EHR_AUDIT_REPO_ROOT}/output/rrt_first.csv
--   ${EHR_AUDIT_REPO_ROOT}/output/rrt_timing_windows_v2.csv
--   ${EHR_AUDIT_REPO_ROOT}/output/rrt_hte_base_times.csv
--   ${EHR_AUDIT_REPO_ROOT}/output/features_0_6.csv
-- =============================================================

CREATE SCHEMA IF NOT EXISTS project_sa_aki;

\echo 'Step 1. Build cohort_icustay table'

DROP TABLE IF EXISTS project_sa_aki.cohort_icustay;
CREATE TABLE project_sa_aki.cohort_icustay AS
SELECT DISTINCT
    f.stay_id,
    i.subject_id,
    i.hadm_id,
    i.intime,
    i.outtime
FROM project_sa_aki.final_full f
INNER JOIN mimiciv_icu.icustays i
    ON f.stay_id = i.stay_id;

CREATE INDEX IF NOT EXISTS idx_cohort_icustay_stay_id
    ON project_sa_aki.cohort_icustay(stay_id);
CREATE INDEX IF NOT EXISTS idx_cohort_icustay_hadm_id
    ON project_sa_aki.cohort_icustay(hadm_id);

\echo 'Step 2. Build first RRT time table'

-- This block dynamically adapts to common mimiciv_derived.rrt column names.
-- If no dialysis indicator column is detected, all rows in mimiciv_derived.rrt
-- with non-missing charttime are treated as RRT events.
DO $$
DECLARE
    conds text[] := ARRAY[]::text[];
    rrt_cond text;
BEGIN
    IF NOT EXISTS (
        SELECT 1
        FROM information_schema.tables
        WHERE table_schema = 'mimiciv_derived'
          AND table_name = 'rrt'
    ) THEN
        RAISE EXCEPTION 'Table mimiciv_derived.rrt not found. Please create derived RRT table first.';
    END IF;

    IF EXISTS (
        SELECT 1 FROM information_schema.columns
        WHERE table_schema = 'mimiciv_derived' AND table_name = 'rrt'
          AND column_name = 'dialysis_active'
    ) THEN
        conds := array_append(conds, $$lower(coalesce(r.dialysis_active::text, '')) IN ('1','t','true','yes','y')$$);
    END IF;

    IF EXISTS (
        SELECT 1 FROM information_schema.columns
        WHERE table_schema = 'mimiciv_derived' AND table_name = 'rrt'
          AND column_name = 'dialysis_present'
    ) THEN
        conds := array_append(conds, $$lower(coalesce(r.dialysis_present::text, '')) IN ('1','t','true','yes','y')$$);
    END IF;

    IF EXISTS (
        SELECT 1 FROM information_schema.columns
        WHERE table_schema = 'mimiciv_derived' AND table_name = 'rrt'
          AND column_name = 'rrt'
    ) THEN
        conds := array_append(conds, $$lower(coalesce(r.rrt::text, '')) IN ('1','t','true','yes','y')$$);
    END IF;

    IF EXISTS (
        SELECT 1 FROM information_schema.columns
        WHERE table_schema = 'mimiciv_derived' AND table_name = 'rrt'
          AND column_name = 'dialysis_type'
    ) THEN
        conds := array_append(conds, $$(r.dialysis_type IS NOT NULL AND r.dialysis_type::text <> '')$$);
    END IF;

    IF array_length(conds, 1) IS NULL THEN
        rrt_cond := 'r.charttime IS NOT NULL';
    ELSE
        rrt_cond := 'r.charttime IS NOT NULL AND (' || array_to_string(conds, ' OR ') || ')';
    END IF;

    EXECUTE 'DROP TABLE IF EXISTS project_sa_aki.rrt_first';

    EXECUTE format($SQL$
        CREATE TABLE project_sa_aki.rrt_first AS
        WITH rrt0 AS (
            SELECT
                c.stay_id,
                r.charttime,
                EXTRACT(EPOCH FROM (r.charttime - c.intime)) / 3600.0 AS rrt_hours
            FROM mimiciv_derived.rrt r
            INNER JOIN project_sa_aki.cohort_icustay c
                ON r.stay_id = c.stay_id
            WHERE %s
        )
        SELECT
            stay_id,
            MIN(rrt_hours) AS first_rrt_hours
        FROM rrt0
        GROUP BY stay_id
        ORDER BY stay_id
    $SQL$, rrt_cond);
END $$;

CREATE INDEX IF NOT EXISTS idx_rrt_first_stay_id
    ON project_sa_aki.rrt_first(stay_id);

\echo 'Step 3. Build RRT timing windows and landmark time variables'

DROP TABLE IF EXISTS project_sa_aki.rrt_hte_base_times;
CREATE TABLE project_sa_aki.rrt_hte_base_times AS
SELECT
    c.stay_id,
    c.subject_id,
    c.hadm_id,
    c.intime,
    c.outtime,
    rf.first_rrt_hours,
    CASE
        WHEN COALESCE(a.deathtime, p.dod::timestamp) IS NULL THEN NULL
        ELSE EXTRACT(EPOCH FROM (COALESCE(a.deathtime, p.dod::timestamp) - c.intime)) / 3600.0
    END AS death_time_hours,
    CASE
        WHEN COALESCE(a.deathtime, p.dod::timestamp) IS NOT NULL
         AND EXTRACT(EPOCH FROM (COALESCE(a.deathtime, p.dod::timestamp) - c.intime)) / 3600.0 <= 24
        THEN 1 ELSE 0
    END AS mortality_24h,
    CASE
        WHEN COALESCE(a.deathtime, p.dod::timestamp) IS NOT NULL
         AND EXTRACT(EPOCH FROM (COALESCE(a.deathtime, p.dod::timestamp) - c.intime)) / 3600.0 <= 72
        THEN 1 ELSE 0
    END AS mortality_72h,
    CASE
        WHEN COALESCE(a.deathtime, p.dod::timestamp) IS NOT NULL
         AND EXTRACT(EPOCH FROM (COALESCE(a.deathtime, p.dod::timestamp) - c.intime)) / 3600.0 <= 24*30
        THEN 1 ELSE 0
    END AS mortality_30d_from_time,
    CASE
        WHEN COALESCE(a.deathtime, p.dod::timestamp) IS NOT NULL
         AND EXTRACT(EPOCH FROM (COALESCE(a.deathtime, p.dod::timestamp) - c.intime)) / 3600.0 <= 24*90
        THEN 1 ELSE 0
    END AS mortality_90d_from_time
FROM project_sa_aki.cohort_icustay c
LEFT JOIN project_sa_aki.rrt_first rf
    ON c.stay_id = rf.stay_id
LEFT JOIN mimiciv_hosp.admissions a
    ON c.hadm_id = a.hadm_id
LEFT JOIN mimiciv_hosp.patients p
    ON c.subject_id = p.subject_id;

CREATE INDEX IF NOT EXISTS idx_rrt_hte_base_times_stay_id
    ON project_sa_aki.rrt_hte_base_times(stay_id);

DROP TABLE IF EXISTS project_sa_aki.rrt_timing_windows_v2;
CREATE TABLE project_sa_aki.rrt_timing_windows_v2 AS
SELECT
    stay_id,
    first_rrt_hours,
    CASE
        WHEN first_rrt_hours IS NULL THEN 'No RRT'
        WHEN first_rrt_hours < 0 THEN 'RRT before ICU'
        WHEN first_rrt_hours >= 0  AND first_rrt_hours < 6  THEN 'RRT 0-6h'
        WHEN first_rrt_hours >= 6  AND first_rrt_hours < 24 THEN 'RRT 6-24h'
        WHEN first_rrt_hours >= 24 AND first_rrt_hours < 72 THEN 'RRT 24-72h'
        WHEN first_rrt_hours >= 72 THEN 'RRT after 72h'
        ELSE 'Unknown'
    END AS rrt_window_v2,
    CASE WHEN first_rrt_hours IS NOT NULL AND first_rrt_hours < 6 THEN 1 ELSE 0 END AS rrt_before_6h,
    CASE WHEN first_rrt_hours IS NOT NULL AND first_rrt_hours >= 6 AND first_rrt_hours < 24 THEN 1 ELSE 0 END AS rrt_6_24h,
    CASE WHEN first_rrt_hours IS NOT NULL AND first_rrt_hours < 24 THEN 1 ELSE 0 END AS rrt_before_24h,
    CASE WHEN first_rrt_hours IS NOT NULL AND first_rrt_hours >= 24 AND first_rrt_hours < 72 THEN 1 ELSE 0 END AS rrt_24_72h,
    CASE WHEN first_rrt_hours IS NOT NULL AND first_rrt_hours >= 72 THEN 1 ELSE 0 END AS rrt_after_72h
FROM project_sa_aki.rrt_hte_base_times;

CREATE INDEX IF NOT EXISTS idx_rrt_timing_windows_v2_stay_id
    ON project_sa_aki.rrt_timing_windows_v2(stay_id);

\echo 'Step 4. Extract 0-6h features for early phenotype'

DROP TABLE IF EXISTS project_sa_aki.features_0_6;
CREATE TABLE project_sa_aki.features_0_6 AS
WITH cohort AS (
    SELECT stay_id, subject_id, hadm_id, intime
    FROM project_sa_aki.cohort_icustay
),
chem AS (
    SELECT
        c.stay_id,
        MAX(ch.creatinine) AS creatinine_max_0_6h,
        MAX(ch.bun) AS bun_max_0_6h,
        MAX(ch.potassium) AS potassium_max_0_6h,
        MIN(ch.bicarbonate) AS bicarbonate_min_0_6h,
        MIN(ch.sodium) AS sodium_min_0_6h,
        MAX(ch.chloride) AS chloride_max_0_6h,
        MAX(ch.glucose) AS glucose_max_0_6h,
        MAX(ch.aniongap) AS aniongap_max_0_6h
    FROM cohort c
    LEFT JOIN mimiciv_derived.chemistry ch
        ON c.hadm_id = ch.hadm_id
       AND ch.charttime >= c.intime
       AND ch.charttime <  c.intime + INTERVAL '6 hours'
    GROUP BY c.stay_id
),
cbc AS (
    SELECT
        c.stay_id,
        MIN(cb.platelet) AS platelets_min_0_6h,
        MAX(cb.wbc) AS wbc_max_0_6h,
        MIN(cb.hemoglobin) AS hemoglobin_min_0_6h
    FROM cohort c
    LEFT JOIN mimiciv_derived.complete_blood_count cb
        ON c.hadm_id = cb.hadm_id
       AND cb.charttime >= c.intime
       AND cb.charttime <  c.intime + INTERVAL '6 hours'
    GROUP BY c.stay_id
),
coag AS (
    SELECT
        c.stay_id,
        MAX(co.inr) AS inr_max_0_6h,
        MAX(co.pt) AS pt_max_0_6h,
        MAX(co.ptt) AS ptt_max_0_6h
    FROM cohort c
    LEFT JOIN mimiciv_derived.coagulation co
        ON c.hadm_id = co.hadm_id
       AND co.charttime >= c.intime
       AND co.charttime <  c.intime + INTERVAL '6 hours'
    GROUP BY c.stay_id
),
enzyme AS (
    SELECT
        c.stay_id,
        MAX(en.bilirubin_total) AS bilirubin_total_max_0_6h,
        MAX(en.alt) AS alt_max_0_6h,
        MAX(en.ast) AS ast_max_0_6h
    FROM cohort c
    LEFT JOIN mimiciv_derived.enzyme en
        ON c.hadm_id = en.hadm_id
       AND en.charttime >= c.intime
       AND en.charttime <  c.intime + INTERVAL '6 hours'
    GROUP BY c.stay_id
),
bg AS (
    SELECT
        c.stay_id,
        MIN(bg.ph) AS ph_min_0_6h,
        MAX(bg.lactate) AS lactate_max_0_6h,
        MIN(bg.po2) AS po2_min_0_6h,
        MAX(bg.pco2) AS pco2_max_0_6h,
        MIN(bg.baseexcess) AS baseexcess_min_0_6h,
        MIN(bg.so2) AS so2_min_0_6h,
        MIN(bg.pao2fio2ratio) AS pao2fio2ratio_min_0_6h
    FROM cohort c
    LEFT JOIN mimiciv_derived.bg bg
        ON c.hadm_id = bg.hadm_id
       AND bg.charttime >= c.intime
       AND bg.charttime <  c.intime + INTERVAL '6 hours'
    GROUP BY c.stay_id
),
vital AS (
    SELECT
        c.stay_id,
        MIN(vs.mbp) AS mbp_min_0_6h,
        MAX(vs.heart_rate) AS heart_rate_max_0_6h,
        MAX(vs.resp_rate) AS resp_rate_max_0_6h,
        MAX(vs.temperature) AS temperature_max_0_6h,
        MIN(vs.spo2) AS spo2_min_0_6h,
        MIN(vs.sbp) AS sbp_min_0_6h
    FROM cohort c
    LEFT JOIN mimiciv_derived.vitalsign vs
        ON c.stay_id = vs.stay_id
       AND vs.charttime >= c.intime
       AND vs.charttime <  c.intime + INTERVAL '6 hours'
    GROUP BY c.stay_id
),
urine AS (
    SELECT
        c.stay_id,
        SUM(uo.urineoutput) AS urine_output_0_6h_ml
    FROM cohort c
    LEFT JOIN mimiciv_derived.urine_output uo
        ON c.stay_id = uo.stay_id
       AND uo.charttime >= c.intime
       AND uo.charttime <  c.intime + INTERVAL '6 hours'
    GROUP BY c.stay_id
),
gcs AS (
    SELECT
        c.stay_id,
        MIN(g.gcs) AS gcs_min_0_6h
    FROM cohort c
    LEFT JOIN mimiciv_derived.gcs g
        ON c.stay_id = g.stay_id
       AND g.charttime >= c.intime
       AND g.charttime <  c.intime + INTERVAL '6 hours'
    GROUP BY c.stay_id
)
SELECT
    c.stay_id,
    chem.creatinine_max_0_6h,
    chem.bun_max_0_6h,
    chem.potassium_max_0_6h,
    chem.bicarbonate_min_0_6h,
    chem.sodium_min_0_6h,
    chem.chloride_max_0_6h,
    chem.glucose_max_0_6h,
    chem.aniongap_max_0_6h,
    cbc.platelets_min_0_6h,
    cbc.wbc_max_0_6h,
    cbc.hemoglobin_min_0_6h,
    coag.inr_max_0_6h,
    coag.pt_max_0_6h,
    coag.ptt_max_0_6h,
    enzyme.bilirubin_total_max_0_6h,
    enzyme.alt_max_0_6h,
    enzyme.ast_max_0_6h,
    bg.ph_min_0_6h,
    bg.lactate_max_0_6h,
    bg.po2_min_0_6h,
    bg.pco2_max_0_6h,
    bg.baseexcess_min_0_6h,
    bg.so2_min_0_6h,
    bg.pao2fio2ratio_min_0_6h,
    vital.mbp_min_0_6h,
    vital.heart_rate_max_0_6h,
    vital.resp_rate_max_0_6h,
    vital.temperature_max_0_6h,
    vital.spo2_min_0_6h,
    vital.sbp_min_0_6h,
    urine.urine_output_0_6h_ml,
    gcs.gcs_min_0_6h
FROM cohort c
LEFT JOIN chem   ON c.stay_id = chem.stay_id
LEFT JOIN cbc    ON c.stay_id = cbc.stay_id
LEFT JOIN coag   ON c.stay_id = coag.stay_id
LEFT JOIN enzyme ON c.stay_id = enzyme.stay_id
LEFT JOIN bg     ON c.stay_id = bg.stay_id
LEFT JOIN vital  ON c.stay_id = vital.stay_id
LEFT JOIN urine  ON c.stay_id = urine.stay_id
LEFT JOIN gcs    ON c.stay_id = gcs.stay_id
ORDER BY c.stay_id;

CREATE INDEX IF NOT EXISTS idx_features_0_6_stay_id
    ON project_sa_aki.features_0_6(stay_id);

\echo 'Step 5. Export CSV files'

\copy (SELECT stay_id, first_rrt_hours FROM project_sa_aki.rrt_first ORDER BY stay_id) TO :'rrt_first_csv' WITH CSV HEADER;

\copy (SELECT * FROM project_sa_aki.rrt_timing_windows_v2 ORDER BY stay_id) TO :'rrt_windows_csv' WITH CSV HEADER;

\copy (SELECT * FROM project_sa_aki.rrt_hte_base_times ORDER BY stay_id) TO :'rrt_base_times_csv' WITH CSV HEADER;

\copy (SELECT * FROM project_sa_aki.features_0_6 ORDER BY stay_id) TO :'features_0_6_csv' WITH CSV HEADER;

\echo 'Step 6. QC summaries'

SELECT 'cohort_icustay' AS table_name, COUNT(*) AS n FROM project_sa_aki.cohort_icustay
UNION ALL
SELECT 'rrt_first', COUNT(*) FROM project_sa_aki.rrt_first
UNION ALL
SELECT 'rrt_hte_base_times', COUNT(*) FROM project_sa_aki.rrt_hte_base_times
UNION ALL
SELECT 'features_0_6', COUNT(*) FROM project_sa_aki.features_0_6;

SELECT
    rrt_window_v2,
    COUNT(*) AS n,
    ROUND(100.0 * COUNT(*) / SUM(COUNT(*)) OVER (), 2) AS pct
FROM project_sa_aki.rrt_timing_windows_v2
GROUP BY rrt_window_v2
ORDER BY
    CASE rrt_window_v2
        WHEN 'RRT before ICU' THEN 1
        WHEN 'RRT 0-6h' THEN 2
        WHEN 'RRT 6-24h' THEN 3
        WHEN 'RRT 24-72h' THEN 4
        WHEN 'RRT after 72h' THEN 5
        WHEN 'No RRT' THEN 6
        ELSE 99
    END;

-- Missingness overview for key 0-6h features.
SELECT
    COUNT(*) AS n,
    ROUND(100.0 * AVG((creatinine_max_0_6h IS NULL)::int), 2) AS miss_creatinine_pct,
    ROUND(100.0 * AVG((bun_max_0_6h IS NULL)::int), 2) AS miss_bun_pct,
    ROUND(100.0 * AVG((potassium_max_0_6h IS NULL)::int), 2) AS miss_potassium_pct,
    ROUND(100.0 * AVG((bicarbonate_min_0_6h IS NULL)::int), 2) AS miss_bicarbonate_pct,
    ROUND(100.0 * AVG((ph_min_0_6h IS NULL)::int), 2) AS miss_ph_pct,
    ROUND(100.0 * AVG((lactate_max_0_6h IS NULL)::int), 2) AS miss_lactate_pct,
    ROUND(100.0 * AVG((urine_output_0_6h_ml IS NULL)::int), 2) AS miss_urine_pct,
    ROUND(100.0 * AVG((mbp_min_0_6h IS NULL)::int), 2) AS miss_mbp_pct,
    ROUND(100.0 * AVG((gcs_min_0_6h IS NULL)::int), 2) AS miss_gcs_pct
FROM project_sa_aki.features_0_6;

\echo 'Done. D4 timing/support input exports completed.'
