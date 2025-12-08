#' @note
#' Needs higher tolerance than macro_fisher_null to get superior LL

fit_null_discrete_pseudohuber <- function(
    Y,
    X,
    k_constr,
    j_constr,
    j_ref,
    trackB = TRUE,
    maxit = 1000, 
    tol = 1e-8,
    ls_max = 20, 
    ls_rho = 0.5, 
    max_step = 1
) {
  
  n <- nrow(Y)
  J <- ncol(Y)
  p <- ncol(X)
  
  ### check we're in the right paradigm
  distinct_xx <- unique(X)
  stopifnot(ncol(X) == nrow(distinct_xx))
  stopifnot(all(X[,1] == 1)) ### need a baseline category to simplify the parametrization
  X_wo_1s <- distinct_xx[, -1]
  
  
  ## Rearrange to have j_ref be the last (J-th) column
  ## Setup: 
  ##    put j_ref to be the J-th column
  ##    put j_constr to be the 1st column
  ## Later, permute cols of B-hat
  Y_orig <- Y
  new_order <- c(j_constr, setdiff(1:J, c(j_constr, j_ref)), j_ref)
  Y <- Y_orig[, new_order]
  
  groups <- split(
    seq_len(nrow(X)),                        # row indices of original data
    apply(X, 1, function(r)
      paste(r, collapse = "_"))               # grouping key by row contents
  )
  
  groups <- groups[order(sapply(groups, min))]
  
  totals <- lapply(groups, function(x) { 
    apply(Y[x, , drop = FALSE], 2, sum)
  })
  
  Y_sum <- do.call(rbind, totals)
  
  ## this should be equivalent to fit_null with j_ref=J
  out <- fit_null_discrete_micro_fs(
    Y = Y_sum, 
    X = X_wo_1s, 
    k_constr=k_constr, 
    j_constr=1, ## because it's been moved to the first column
    constraint_fn=function(x) {  pseudohuber_median(c(x, 0)) },       
    constraint_grad_fn= function(x) {  x <- radEmu::dpseudohuber_median_dx(c(x, 0)); x[-length(x)]},
    maxit = maxit, 
    tol = tol,
    ls_max = ls_max, 
    ls_rho = ls_rho, 
    max_step = max_step
  )

  B <- out$B
  
  z <- update_z(Y, X, B)
  log_means <- X %*% B + matrix(z, ncol = 1) %*% matrix(1, ncol = J, nrow = 1)
  ll_new <- sum(Y * log_means - exp(log_means))
  
  ### now reorder B to match what was input 
  return(list(
    "B" = B[, order(new_order)],
    "k_constr" = k_constr,
    "j_constr" = j_constr,
    "niter" = out$iter,
    "loglik" = ll_new,
    "it_df" = out$history,
    "converged" = (out$iter < maxit)
  ))
  
}