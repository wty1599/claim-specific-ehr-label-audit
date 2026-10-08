-- Read-only dictionary screen. Run against the authorized local MIMIC-IV database.
-- This query does not inspect patient records or produce treatment/outcome counts.
-- Its output is a CANDIDATE list, not an approved IV itemid list: review each
-- category and the associated inputevents order fields before freezing IDs.

SELECT
    itemid,
    label,
    abbreviation,
    category,
    unitname,
    param_type,
    linksto
FROM mimiciv_icu.d_items
WHERE linksto = 'inputevents'
  AND (
      label ~* '(furosemide|lasix|bumetanide|bumex|torsemide|demadex)'
      OR COALESCE(abbreviation, '') ~* '(furosemide|lasix|bumetanide|bumex|torsemide|demadex)'
  )
ORDER BY LOWER(label), itemid;
