-- ============================================================



-- ============================================================
\set ON_ERROR_STOP on
\pset pager off
\if :{?eicu_external_csv}
\else
\set eicu_external_csv './outputs/sql/eicu_external.csv'
\endif


DROP TABLE IF EXISTS public.saaki_cohort_dedup;
CREATE TABLE public.saaki_cohort_dedup AS
SELECT DISTINCT ON (patientunitstayid) *
FROM public.saaki_cohort
ORDER BY patientunitstayid, apache_iva DESC NULLS LAST;

\echo '===== Deduplicated cohort ====='
SELECT COUNT(*) total, COUNT(DISTINCT patientunitstayid) uniq FROM public.saaki_cohort_dedup;


DROP TABLE public.saaki_cohort;
ALTER TABLE public.saaki_cohort_dedup RENAME TO saaki_cohort;

\echo '===== AKI stage distribution (expected unchanged) ====='
SELECT aki_stage_0_24h, COUNT(*) FROM public.saaki_cohort GROUP BY 1 ORDER BY 1;
\echo '===== In-hospital mortality (historical expectation: unchanged ~0.236) ====='
SELECT ROUND(AVG(hosp_mortality)::numeric,3) FROM public.saaki_cohort;


DROP TABLE IF EXISTS public.eicu_external_dedup;
CREATE TABLE public.eicu_external_dedup AS
SELECT DISTINCT ON (patientunitstayid) *
FROM public.eicu_external
ORDER BY patientunitstayid, apache_iva DESC NULLS LAST;

DROP TABLE public.eicu_external;
ALTER TABLE public.eicu_external_dedup RENAME TO eicu_external;

\echo '===== After deduplication eicu_external ====='
SELECT COUNT(*) total, COUNT(DISTINCT patientunitstayid) uniq FROM public.eicu_external;


\copy (SELECT * FROM public.eicu_external) TO :'eicu_external_csv' CSV HEADER
\echo 'Re-exported eicu_external.csv (After deduplication)'
