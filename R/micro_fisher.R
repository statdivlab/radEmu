#fisher scoring update for j-th column of B in unconstrained optimization
#z held constant
micro_fisher <- function(X, Yj, Bj, z, stepsize = 1, c1 = 0.1) {
  log_means <- X %*% Bj + z
  means <- exp(log_means)

  #info in Bj
  info <- t(X) %*% diag(as.numeric(means)) %*% X
  lj_grad <- colSums(diag(as.numeric(Yj - means)) %*% X)

  #make update a try-error to start
  update <- try(stop(), silent = TRUE)

  if (nrow(info) > 1) {
    info_avg_diag <- diag(rep(sqrt(mean(diag(info)^2)), nrow(info)))
  } else {
    info_avg_diag <- abs(info)
  }

  #try to compute update direction as is, but if we run into numerical
  #invertibility issues, regularize info and try again
  regularization <- 0
  while (inherits(update, "try-error")) {
    update <- try(
      qr.solve(info + regularization * info_avg_diag, lj_grad),
      silent = TRUE
    )
    regularization <- ifelse(regularization == 0, 0.01, 10 * regularization)
  }

  #use armijo rule to choose step size
  obj <- -sum(Yj * log_means - means)
  obj_grad <- -lj_grad

  suff_decrease_term <- c1 * sum(obj_grad * update)

  suff_decrease <- FALSE
  while (!(suff_decrease)) {
    prop_Bj <- Bj + stepsize * update
    prop_log_mu <- X %*% prop_Bj + z
    prop_obj <- -sum(Yj * prop_log_mu - exp(prop_log_mu))

    suff_decrease <- prop_obj <= obj + suff_decrease_term * stepsize

    stepsize <- stepsize * 0.5
  }

  stepsize <- stepsize / 0.5

  return(stepsize * update)
}
