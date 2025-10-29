null_repar_ll <- function(
  x,
  js,
  B,
  z,
  p,
  Y,
  X,
  j_constr,
  k_constr,
  constraint_fn,
  ref_set = NULL
) {
  Bjs <- B[, c(js, j_constr), drop = FALSE]
  njs <- length(js)
  for (jind in 1:njs) {
    Bjs[, jind] <- Bjs[, jind] + x[1:p + (jind - 1) * p]
  }
  Bjs[-k_constr, njs + 1] <- Bjs[-k_constr, njs + 1] + x[p * njs + 1:(p - 1)]
  Bk_temp <- B[k_constr, ]
  Bk_temp[js] <- Bjs[k_constr, 1:njs]
  Bk_temp[j_constr] <- NA
  Bjs[k_constr, njs + 1] <- compute_constraint_value(constraint_fn, Bk_temp, j_constr, ref_set)
  #
  log_means_js <- X %*% Bjs
  #
  for (i in 1:nrow(Y)) {
    log_means_js[i, ] <- log_means_js[i, ] + z[i]
  }
  #
  return(sum(Y[, c(js, j_constr)] * log_means_js - exp(log_means_js)))
}
