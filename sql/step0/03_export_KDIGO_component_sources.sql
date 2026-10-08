\set ON_ERROR_STOP on

-- Descriptive re-export only.
-- This query uses the same ICU-admission 0-24 h window and component sources
-- as sql/build_baseline_covars_v2.sql. It does not redefine KDIGO.

CREATE TEMP TABLE kdigo_component_sources_0_24h_audit AS
WITH aki_24h AS (
  SELECT
    ie.stay_id,
    MAX(k.aki_stage_creat) AS stage_creat_0_24h,
    MAX(k.aki_stage_uo) AS stage_uo_0_24h
  FROM mimiciv_icu.icustays AS ie
  LEFT JOIN mimiciv_derived.kdigo_stages AS k
    ON k.stay_id = ie.stay_id
   AND k.charttime >= ie.intime
   AND k.charttime < ie.intime + INTERVAL '24 hour'
  GROUP BY ie.stay_id
),
rrt_24h AS (
  SELECT
    ie.stay_id,
    MAX(
      CASE
        WHEN r.dialysis_active = 1 THEN 1
        ELSE 0
      END
    ) AS rrt_0_24h
  FROM mimiciv_icu.icustays AS ie
  LEFT JOIN mimiciv_derived.rrt AS r
    ON r.stay_id = ie.stay_id
   AND r.charttime >= ie.intime
   AND r.charttime < ie.intime + INTERVAL '24 hour'
  GROUP BY ie.stay_id
),
component_export AS (
  SELECT
    b.stay_id,
    COALESCE(a.stage_creat_0_24h, 0)::integer AS stage_creat_0_24h,
    COALESCE(a.stage_uo_0_24h, 0)::integer AS stage_uo_0_24h,
    COALESCE(r.rrt_0_24h, 0)::integer AS rrt_0_24h,
    b.aki_stage_0_24h::integer AS aki_stage_0_24h_existing,
    GREATEST(
      COALESCE(a.stage_creat_0_24h, 0),
      COALESCE(a.stage_uo_0_24h, 0),
      CASE WHEN COALESCE(r.rrt_0_24h, 0) = 1 THEN 3 ELSE 0 END
    )::integer AS aki_stage_0_24h_reconstructed
  FROM project_sa_aki.baseline_covars AS b
  LEFT JOIN aki_24h AS a
    ON a.stay_id = b.stay_id
  LEFT JOIN rrt_24h AS r
    ON r.stay_id = b.stay_id
)
SELECT
  stay_id,
  stage_creat_0_24h,
  stage_uo_0_24h,
  rrt_0_24h,
  aki_stage_0_24h_existing,
  aki_stage_0_24h_reconstructed,
  (aki_stage_0_24h_existing = aki_stage_0_24h_reconstructed)::integer
    AS stage_reconstruction_match
FROM component_export
ORDER BY stay_id;

\copy (SELECT * FROM kdigo_component_sources_0_24h_audit ORDER BY stay_id) TO '/path/to/MIMIC-IV/kdigo_criterion_decomposition_20260725/data_extract/kdigo_component_sources_0_24h.csv' WITH (FORMAT CSV, HEADER TRUE);

\echo 'KDIGO component-source export completed.'
\echo 'No existing cohort, KDIGO stage, MICE object, or cluster label was modified.'
