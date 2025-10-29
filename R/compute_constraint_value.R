compute_constraint_value <- function(constraint_fn, vec, j_constr, ref_set = NULL) {
  if (is.null(ref_set)) {
    return(constraint_fn(vec[-j_constr]))
  }
  if (attr(constraint_fn, "constraint_type") == "symmetric_subset:pseudohuber") {
    return(pseudohuber_median(vec[setdiff(ref_set, j_constr)]))
  } else if (attr(constraint_fn, "constraint_type") == "symmetric_subset:mean") {
    return(mean(vec[setdiff(ref_set, j_constr)]))
  } else {
    return(constraint_fn(vec[-j_constr]))
  }
}