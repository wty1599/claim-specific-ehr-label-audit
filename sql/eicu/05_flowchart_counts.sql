-- ============================================================

-- ============================================================
\pset pager off
\echo '===== eICU cohort flow counts ====='
SELECT 'A. All unit stays'        AS step, COUNT(*) n FROM patient
UNION ALL SELECT 'B. Adults(>=18)', COUNT(*) FROM patient WHERE age='> 89' OR (age ~ '^[0-9]+$' AND age::int>=18)
UNION ALL SELECT 'C. First ICU stay(unitvisitnumber=1)', COUNT(*) FROM patient WHERE unitvisitnumber=1 AND (age='> 89' OR (age ~ '^[0-9]+$' AND age::int>=18))
UNION ALL SELECT 'D. Sepsis (diagnosis text)', COUNT(DISTINCT patientunitstayid) FROM diagnosis WHERE lower(diagnosisstring) LIKE '%sepsis%' OR lower(diagnosisstring) LIKE '%septic%'
UNION ALL SELECT 'E. Final sepsis cohort after prior CKD/ESRD exclusions', COUNT(*) FROM public.saaki_cohort;
