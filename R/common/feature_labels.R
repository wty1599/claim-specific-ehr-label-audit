## ============================================================




##   FEAT_UNIT["ph_min"]         -> ""  ; FEAT_SYSTEM["ph_min"] -> "Acid–base"
## ============================================================

FEAT_LABEL <- c(
  mbp_min             = "Mean arterial pressure (min)",
  heart_rate_max      = "Heart rate (max)",
  resp_rate_max       = "Respiratory rate (max)",
  temperature_max     = "Temperature (max)",
  temperature_min     = "Temperature (min)",
  lactate_max         = "Lactate (max)",
  glucose_max         = "Glucose (max)",
  ph_min              = "Arterial pH (min)",
  pco2_max            = "PaCO2 (max)",
  bicarbonate_min     = "Bicarbonate (min)",
  aniongap_max        = "Anion gap (max)",
  pao2fio2ratio_min   = "PaO2/FiO2 ratio (min)",
  spo2_min            = "SpO2 (min)",
  creatinine_max      = "Serum creatinine (max)",
  bun_max             = "Blood urea nitrogen (max)",
  urine_output_24h_ml = "24-hour urine output",
  gcs_min             = "Glasgow Coma Scale (min)",
  bilirubin_total_max = "Total bilirubin (max)",
  alt_max             = "ALT (max)",
  inr_max             = "INR (max)",
  ptt_max             = "aPTT (max)",
  wbc_max             = "White blood cell count (max)",
  platelets_min       = "Platelet count (min)",
  hemoglobin_min      = "Hemoglobin (min)",
  abs_lymphocytes_min = "Absolute lymphocyte count (min)",
  sodium_min          = "Sodium (min)",
  sodium_max          = "Sodium (max)",
  potassium_min       = "Potassium (min)",
  potassium_max       = "Potassium (max)",
  chloride_min        = "Chloride (min)",
  chloride_max        = "Chloride (max)",
  calcium_min         = "Calcium (min)",
  calcium_max         = "Calcium (max)"
)

FEAT_UNIT <- c(
  mbp_min="mmHg", heart_rate_max="beats/min", resp_rate_max="breaths/min",
  temperature_max="\u00B0C", temperature_min="\u00B0C", lactate_max="mmol/L",
  glucose_max="mg/dL", ph_min="", pco2_max="mmHg", bicarbonate_min="mmol/L",
  aniongap_max="mmol/L", pao2fio2ratio_min="mmHg", spo2_min="%",
  creatinine_max="mg/dL", bun_max="mg/dL", urine_output_24h_ml="mL",
  gcs_min="points", bilirubin_total_max="mg/dL", alt_max="U/L", inr_max="",
  ptt_max="s", wbc_max="\u00D710\u2079/L", platelets_min="\u00D710\u2079/L",
  hemoglobin_min="g/dL", abs_lymphocytes_min="\u00D710\u2079/L",
  sodium_min="mmol/L", sodium_max="mmol/L", potassium_min="mmol/L",
  potassium_max="mmol/L", chloride_min="mmol/L", chloride_max="mmol/L",
  calcium_min="mg/dL", calcium_max="mg/dL"
)

FEAT_SYSTEM <- c(
  mbp_min="Circulation/Vitals", heart_rate_max="Circulation/Vitals",
  resp_rate_max="Circulation/Vitals", temperature_max="Circulation/Vitals",
  temperature_min="Circulation/Vitals", lactate_max="Perfusion/Metabolic",
  glucose_max="Perfusion/Metabolic", ph_min="Acid-base", pco2_max="Acid-base",
  bicarbonate_min="Acid-base", aniongap_max="Acid-base",
  pao2fio2ratio_min="Respiratory", spo2_min="Respiratory",
  creatinine_max="Renal", bun_max="Renal", urine_output_24h_ml="Renal",
  gcs_min="Neurologic", bilirubin_total_max="Hepatic", alt_max="Hepatic",
  inr_max="Coagulation", ptt_max="Coagulation", wbc_max="Hematologic",
  platelets_min="Hematologic", hemoglobin_min="Hematologic",
  abs_lymphocytes_min="Hematologic", sodium_min="Electrolytes",
  sodium_max="Electrolytes", potassium_min="Electrolytes",
  potassium_max="Electrolytes", chloride_min="Electrolytes",
  chloride_max="Electrolytes", calcium_min="Electrolytes", calcium_max="Electrolytes"
)


FEAT_LABEL_FULL <- vapply(names(FEAT_LABEL), function(k){
  u <- FEAT_UNIT[[k]]; if(is.null(u)||u=="") FEAT_LABEL[[k]] else paste0(FEAT_LABEL[[k]], ", ", u)
}, character(1)); names(FEAT_LABEL_FULL) <- names(FEAT_LABEL)


lab_feat      <- function(v) ifelse(v %in% names(FEAT_LABEL),      FEAT_LABEL[v],      v)
lab_feat_full <- function(v) ifelse(v %in% names(FEAT_LABEL_FULL), FEAT_LABEL_FULL[v], v)
