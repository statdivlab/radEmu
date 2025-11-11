constraint_grad_vec <- function(constraint_fn,
                                constraint_grad_fn,
                                js_used,
                                Bk_constr,
                                j_constr,
                                p,
                                ref_set = NULL){
  
  if (is.null(ref_set)) {
    # gradient over all except the constrained element
    indices <- setdiff(seq_along(Bk_constr), j_constr)
    grad_full <- constraint_grad_fn(Bk_constr[indices])
  } else {
    indices <- setdiff(ref_set, j_constr)
    if (attr(constraint_fn, "constraint_type") == "symmetric_subset:pseudohuber") {
      grad_full <- dpseudohuber_median_dx(Bk_constr[indices])
    } else if (attr(constraint_fn, "constraint_type") == "symmetric_subset:mean"){
      grad_full <- rep(1/length(indices), length(indices))
    } else {
      grad_full <- constraint_grad_fn(Bk_constr[indices])
    }
  }
  
  # map to js_used positions
  grad_restricted <- sapply(js_used, function(j) {
    idx <- which(indices == j)
    if (length(idx) == 0) return(0)
    grad_full[idx]
  })
  
  return(grad_restricted)
}