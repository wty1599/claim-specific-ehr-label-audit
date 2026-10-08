## Strict schema and decision layer for the governed Domain 1 Lloyd chain.
## This wrapper validates one already-computed repeat. It does not calculate
## scientific diagnostics and is not itself a scientific-result authority.

D1_DECISION_RULE_VERSION <- "D1_LLOYD_THREE_STATE_V2_20260720"
D1_REQUIRED_SHAPE_SOURCE <- "heldout_fisher_discriminant_bh6"
D1_REQUIRED_REFERENCE_SCOPE <- "processed_space_fixed_k2_n_specific"
D1_REQUIRED_PARTITION_SOURCE <- paste(
  "Lloyd K2 fitted on training data; nearest-centroid assignment on",
  "held-out evaluation"
)

d1_not_evaluated <- function(x, code, detail) {
  list(
    state = "NOT_EVALUATED",
    components = x,
    separation_alert = NA,
    shape_alert = NA,
    failure = code,
    failure_code = code,
    failure_detail = detail,
    decision_rule_version = D1_DECISION_RULE_VERSION
  )
}

d1_scalar_finite_numeric <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
}

d1_scalar_nonmissing_logical <- function(x) {
  is.logical(x) && length(x) == 1L && !is.na(x)
}

d1_scalar_nonmissing_character <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

audit_discreteness_v2 <- function(x = NULL) {
  if (is.null(x)) {
    return(d1_not_evaluated(x, "D1_INPUT_NULL", "Input object is NULL."))
  }
  if (!is.list(x)) {
    return(d1_not_evaluated(
      x, "D1_INPUT_NOT_LIST", "Input must be a named list or one-row data frame."
    ))
  }

  required <- c(
    "separation", "gate", "shape_reject", "shape_source",
    "reference_scope", "partition_source", "kmeans_algorithm",
    "kmeans_nstart", "kmeans_itermax"
  )
  if (anyDuplicated(names(x))) {
    return(d1_not_evaluated(
      x,
      "D1_DUPLICATE_COMPONENT",
      "Input contains duplicated component names."
    ))
  }
  missing_names <- setdiff(required, names(x))
  if (length(missing_names)) {
    return(d1_not_evaluated(
      x,
      "D1_COMPONENT_MISSING",
      paste0("Missing required component(s): ", paste(missing_names, collapse = ", "), ".")
    ))
  }

  if (!d1_scalar_finite_numeric(x$separation)) {
    return(d1_not_evaluated(
      x, "D1_SEPARATION_INVALID", "separation must be one finite numeric value."
    ))
  }
  if (!d1_scalar_finite_numeric(x$gate)) {
    return(d1_not_evaluated(
      x, "D1_GATE_INVALID", "gate must be one finite numeric value."
    ))
  }
  if (!d1_scalar_nonmissing_logical(x$shape_reject)) {
    return(d1_not_evaluated(
      x, "D1_SHAPE_REJECT_INVALID",
      "shape_reject must be one non-missing logical value."
    ))
  }
  if (!d1_scalar_nonmissing_character(x$shape_source) ||
      x$shape_source != D1_REQUIRED_SHAPE_SOURCE) {
    return(d1_not_evaluated(
      x, "D1_SHAPE_SOURCE_MISMATCH",
      paste0("shape_source must equal '", D1_REQUIRED_SHAPE_SOURCE, "'.")
    ))
  }
  if (!d1_scalar_nonmissing_character(x$reference_scope) ||
      x$reference_scope != D1_REQUIRED_REFERENCE_SCOPE) {
    return(d1_not_evaluated(
      x, "D1_REFERENCE_SCOPE_MISMATCH",
      paste0("reference_scope must equal '", D1_REQUIRED_REFERENCE_SCOPE, "'.")
    ))
  }
  if (!d1_scalar_nonmissing_character(x$partition_source) ||
      x$partition_source != D1_REQUIRED_PARTITION_SOURCE) {
    return(d1_not_evaluated(
      x, "D1_PARTITION_SOURCE_MISMATCH",
      paste0("partition_source must equal '", D1_REQUIRED_PARTITION_SOURCE, "'.")
    ))
  }
  if (!d1_scalar_nonmissing_character(x$kmeans_algorithm) ||
      x$kmeans_algorithm != "Lloyd") {
    return(d1_not_evaluated(
      x, "D1_KMEANS_ALGORITHM_MISMATCH",
      "kmeans_algorithm must equal 'Lloyd'."
    ))
  }
  if (!d1_scalar_finite_numeric(x$kmeans_nstart) ||
      x$kmeans_nstart != 25) {
    return(d1_not_evaluated(
      x, "D1_KMEANS_NSTART_MISMATCH", "kmeans_nstart must equal 25."
    ))
  }
  if (!d1_scalar_finite_numeric(x$kmeans_itermax) ||
      x$kmeans_itermax != 100) {
    return(d1_not_evaluated(
      x, "D1_KMEANS_ITERMAX_MISMATCH", "kmeans_itermax must equal 100."
    ))
  }

  separation_alert <- x$separation >= x$gate
  shape_alert <- x$shape_reject
  state <- if (separation_alert && shape_alert) {
    "DISCRETE_EVIDENCE"
  } else if (!separation_alert && !shape_alert) {
    "NO_DISCRETE_EVIDENCE"
  } else {
    "INCONCLUSIVE"
  }

  components <- x
  components$separation_alert <- separation_alert
  components$shape_alert <- shape_alert
  list(
    state = state,
    components = components,
    separation_alert = separation_alert,
    shape_alert = shape_alert,
    failure = NULL,
    failure_code = NULL,
    failure_detail = NULL,
    decision_rule_version = D1_DECISION_RULE_VERSION
  )
}
