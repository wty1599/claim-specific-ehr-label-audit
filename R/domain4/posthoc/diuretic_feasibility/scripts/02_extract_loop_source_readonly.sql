-- Read-only local PostgreSQL extraction. Invoke psql with -qAt -F '|'
-- and -o private_not_for_release/loop_source_by_stay_private.psv.
-- The query covers the locked study cohort; the original D4 RDS provides
-- the authoritative hour-24 risk set after an exact stay-id join.

WITH cohort AS (
  SELECT DISTINCT f.stay_id, i.subject_id, i.hadm_id, i.intime, i.outtime
  FROM project_sa_aki.final_full AS f
  JOIN mimiciv_icu.icustays AS i ON i.stay_id = f.stay_id
),
chf AS (
  SELECT subject_id, hadm_id,
         MAX(congestive_heart_failure) AS chf_code,
         COUNT(*) AS charlson_rows
  FROM mimiciv_derived.charlson
  GROUP BY subject_id, hadm_id
),
loop_rows AS (
  SELECT c.stay_id, c.intime, c.outtime,
         e.itemid, e.starttime, e.endtime
  FROM cohort AS c
  JOIN mimiciv_icu.inputevents AS e ON e.stay_id = c.stay_id
  WHERE e.itemid IN (221794, 228340, 229639)
    AND e.amount > 0
    AND e.ordercategoryname IN ('01-Drips', '05-Med Bolus')
    AND e.statusdescription IN
        ('Bolus', 'ChangeDose/Rate', 'FinishedRunning', 'Paused', 'Stopped')
    AND e.starttime IS NOT NULL
    AND e.endtime IS NOT NULL
    AND e.endtime >= e.starttime
    AND e.endtime > c.intime
    AND e.starttime < c.outtime
),
doses AS (
  SELECT stay_id,
         MIN(starttime) AS first_qualifying_start,
         MIN(starttime) FILTER (
           WHERE starttime >= intime + interval '24 hours'
             AND starttime < intime + interval '72 hours'
         ) AS first_start_24_72,
         COUNT(*) AS qualifying_rows,
         COUNT(*) FILTER (
           WHERE starttime < intime + interval '24 hours'
             AND endtime > intime
         ) AS pre24_rows,
         COUNT(*) FILTER (
           WHERE starttime < intime AND endtime > intime
         ) AS preicu_overlap_rows,
         COUNT(*) FILTER (
           WHERE starttime >= intime + interval '24 hours'
             AND starttime < intime + interval '72 hours'
         ) AS starts_24_72_rows
  FROM loop_rows
  GROUP BY stay_id
)
SELECT c.stay_id, c.intime, c.outtime,
       d.first_qualifying_start, d.first_start_24_72,
       COALESCE(d.qualifying_rows, 0),
       COALESCE(d.pre24_rows, 0),
       COALESCE(d.preicu_overlap_rows, 0),
       COALESCE(d.starts_24_72_rows, 0),
       chf.chf_code,
       chf.charlson_rows
FROM cohort AS c
LEFT JOIN doses AS d ON d.stay_id = c.stay_id
LEFT JOIN chf ON chf.subject_id = c.subject_id AND chf.hadm_id = c.hadm_id
ORDER BY c.stay_id;
