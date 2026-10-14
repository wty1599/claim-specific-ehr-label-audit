-- ============================================================
\if :{?final_full_csv}
\else
\set final_full_csv './outputs/sql/final_full.csv'
\endif



-- ============================================================
SET search_path TO project_sa_aki;

DROP TABLE IF EXISTS project_sa_aki.final_full_new;

CREATE TABLE project_sa_aki.final_full_new AS
SELECT
    f.*,
    m.rrt_30d,
    m.persistent_rd_30d,
    m.terminal_scr,
    m.terminal_scr_available,
    m.terminal_scr_source,
    m.make30
FROM project_sa_aki.final_full f
LEFT JOIN project_sa_aki.make30 m ON m.stay_id = f.stay_id;


DROP TABLE project_sa_aki.final_full;
ALTER TABLE project_sa_aki.final_full_new RENAME TO final_full;
CREATE INDEX IF NOT EXISTS idx_final_full_stay ON project_sa_aki.final_full (stay_id);


SELECT
    (SELECT COUNT(*) FROM project_sa_aki.final_full)  AS n_rows,
    (SELECT COUNT(*) FROM information_schema.columns
      WHERE table_schema='project_sa_aki' AND table_name='final_full') AS n_cols,
    (SELECT SUM(make30) FROM project_sa_aki.final_full) AS n_make30;

\copy (SELECT * FROM project_sa_aki.final_full) TO :'final_full_csv' WITH CSV HEADER
