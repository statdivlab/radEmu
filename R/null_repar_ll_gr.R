null_repar_ll_gr <- function(
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
  constraint_grad_fn,
  return_hess = FALSE,
  return_info_inv = TRUE,
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

  Yjs <- Y[, c(js, j_constr)]
  #
  gr <- lapply(1:(njs + 1), function(jind) {
    Matrix::crossprod(
      X,
      Yjs[, jind, drop = FALSE] -
        exp(log_means_js[, jind, drop = FALSE])
    )
  })
  cg <- constraint_grad_vec(
    constraint_grad_fn,
    js_used = js,
    Bk_constr = B[k_constr, ],
    j_constr = j_constr,
    p = ncol(X),
    ref_set = ref_set
  )

  cg_gr_multiplier <- gr[[length(gr)]][k_constr]

  for (jind in 1:njs) {
    gr[[jind]][k_constr] <- gr[[jind]][k_constr] +
      cg[jind] * cg_gr_multiplier
  }
  gr[[length(gr)]] <- gr[[length(gr)]][-k_constr]
  gr <- do.call(c, gr)

  if (return_hess) {
    info_diags <- lapply(
      1:njs,
      function(jind) Matrix::crossprod(X, diag(exp(log_means_js[, jind]))) %*% X
    )

    # old code that was causing a bug with incorrect dimension of info for p > 2
    # info_diags <- c(
    #   info_diags,
    #   Matrix::crossprod(
    #     X[, -k_constr, drop = FALSE],
    #     diag(exp(log_means_js[, njs + 1]))
    #   ) %*%
    #     X[, -k_constr, drop = FALSE]
    # )
    
    # new code to avoid this bug
    curr_len <- length(info_diags)
    info_diags[[curr_len + 1]] <- Matrix::crossprod(
      X[, -k_constr, drop = FALSE], diag(exp(log_means_js[, njs + 1]))) %*%
      X[, -k_constr, drop = FALSE]
    
    info <- Matrix::bdiag(info_diags)
    
    return_obj <- list(
      gr = gr,
      info = info,
      cg = cg,
      cg_info_mult = Matrix::crossprod(
        X[, k_constr, drop = FALSE],
        diag(exp(log_means_js[, njs + 1]))
      ) %*%
        X[, k_constr, drop = FALSE]
    )
    if (return_info_inv) {
      info_inv <- Matrix::bdiag(lapply(info_diags, qr.solve))
      return_obj$info_inv <- info_inv
    }

    return(return_obj)
  } else {
    return(gr)
  }
}
