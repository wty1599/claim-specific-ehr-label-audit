-- ============================================================
-- project_sa_aki : adult first-ICU source cohort
-- ============================================================
-- This is the exact cohort-construction block recovered from the
-- production script project_sa_aki_extract_base.sql (2026-06-29).
-- Provenance limitation: the original standalone producer was not present
-- in the audited source roots; this block was recovered from the prior
-- submission-code release and therefore requires author verification.
-- It creates one first ICU stay per adult subject and is the upstream
-- source for the submitted strict Sepsis-3 cohort.
-- ============================================================

CREATE SCHEMA IF NOT EXISTS project_sa_aki;

SET search_path TO project_sa_aki, mimiciv_hosp, mimiciv_icu, public;

DROP TABLE IF EXISTS project_sa_aki.adult_icu_first CASCADE;

CREATE TABLE project_sa_aki.adult_icu_first AS
WITH icu0 AS (
    SELECT
        p.subject_id,
        p.gender,
        p.anchor_age AS age,
        p.anchor_year,
        p.dod,
        a.hadm_id,
        a.admittime,
        a.dischtime,
        a.deathtime,
        a.hospital_expire_flag,
        a.admission_type,
        a.admission_location,
        a.discharge_location,
        a.insurance,
        a.language,
        a.marital_status,
        a.race,
        i.stay_id,
        i.first_careunit,
        i.last_careunit,
        i.intime,
        i.outtime,
        i.los,
        ROW_NUMBER() OVER (
            PARTITION BY p.subject_id
            ORDER BY i.intime
        ) AS icu_seq
    FROM mimiciv_hosp.patients p
    INNER JOIN mimiciv_hosp.admissions a
        ON p.subject_id = a.subject_id
    INNER JOIN mimiciv_icu.icustays i
        ON a.hadm_id = i.hadm_id
    WHERE p.anchor_age >= 18
)
SELECT *
FROM icu0
WHERE icu_seq = 1;

CREATE INDEX IF NOT EXISTS idx_adult_icu_first_subject
    ON project_sa_aki.adult_icu_first(subject_id);

CREATE INDEX IF NOT EXISTS idx_adult_icu_first_hadm
    ON project_sa_aki.adult_icu_first(hadm_id);

CREATE INDEX IF NOT EXISTS idx_adult_icu_first_stay
    ON project_sa_aki.adult_icu_first(stay_id);

ANALYZE project_sa_aki.adult_icu_first;

SELECT
    COUNT(*) AS n_rows,
    COUNT(DISTINCT subject_id) AS n_subjects,
    COUNT(DISTINCT stay_id) AS n_stays
FROM project_sa_aki.adult_icu_first;
