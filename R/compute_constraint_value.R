compute_constraint_value <- function(constraint_fn, vec, j_constr, ref_set = NULL) {
  if (is.null(ref_set)) {
    return(constraint_fn(vec[-j_constr]))
  }
  if (any(grepl("pseudohuber_median", deparse(body(constraint_fn)), fixed = TRUE))) {
    return(pseudohuber_median(vec[setdiff(ref_set, j_constr)]))
  } else {
    return(mean(vec[setdiff(ref_set, j_constr)]))
  }
}