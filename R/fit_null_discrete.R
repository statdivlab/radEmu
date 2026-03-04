#' fits model with B_kj constrained to equal g(B_k) for constraint fn g, for a symmetric constraint with a discrete design
#' 
#' @param Y Y (with augmentations)
#' @param X design matrix
#' @param k_constr row index of B to constrain
#' @param j_constr col index of B to constrain
#' @param j_ref column index of convenience constraint
#' @param trackB track value of beta across iterations and return? Default is `FALSE`.
#' @param maxit maximum iterations. Default is `5000`.
#' @param tol tolerance for stopping criteria. Algorithm stops when the root mean of the
#' norm of the score vector is less than the tolerance. Default is `0.01`.
#' @param verbose should the algorithm print updates for you? Default is `FALSE`.
#' @param ls_max maximum number of iterations in  the line search. Default is `20`. 
#' @param ls_rho scaling factor in the line search. Default is `0.5`.
#' @param max_step step capping after the line sesarch. Default is `1`. 
#' @param constraint What type of symmetric constraint do we have? Options are `"mean"` and `"pseudohuber`.
#' 
#' @return A list containing elements `B`, `k_constr`, `j_constr`, `niter`
#' and `Bs`. `B` is a matrix containing parameter estimates
#' under the null (obtained by maximum likelihood on augmented observations Y),
#' `k_constr`, and `j_constr` give row and column indexes of the parameter
#' fixed to be equal to the constraint function \eqn{g()} under the null. `niter` is a
#' scalar giving total number of outer iterations used to fit the null model, and
#' `Bs` is a data frame containing values of B by iteration if `trackB` was set
#' equal to TRUE (otherwise it contains a NULL value). 
#' 
fit_null_discrete <- function(
    Y,
    X,
    k_constr,
    j_constr,
    j_ref,
    trackB = FALSE,
    maxit = 5000, 
    tol = 0.1,
    verbose = FALSE,
    ls_max = 20, 
    ls_rho = 0.5, 
    max_step = 1,
    constraint = "pseudohuber"
) {
  
  n <- nrow(Y)
  J <- ncol(Y)
  p <- ncol(X)
  
  ### check we're in the right paradigm
  distinct_xx <- unique(X)
  stopifnot(ncol(X) == nrow(distinct_xx))
  stopifnot(all(X[,1] == 1)) ### need a baseline category to simplify the parametrization
  X_wo_1s <- distinct_xx[, -1]
  
  # reassign j_ref if it is equal to j_constr
  arch_j_ref <- NULL
  if (j_ref == j_constr) {
    arch_j_ref <- j_ref
    # pick first column not equal to j_constr
    j_ref <- setdiff(seq_len(ncol(Y)), j_constr)[1]
  }
  
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
  
  if (constraint == "pseudohuber") {
    constraint_fn <- function(x) {  pseudohuber_median(c(x, 0)) }
    constraint_grad_fn <- function(x) {  x <- radEmu::dpseudohuber_median_dx(c(x, 0)); x[-length(x)]}
  } else if (constraint == "mean") {
    constraint_fn <- function(x) { mean(c(x, 0)) }
    constraint_grad_fn <- function(x) {  x <- rep(1 / (length(x) + 1), length(x))}
  }
  
  ## this should be equivalent to fit_null with j_ref=J
  out <- fit_null_discrete_micro_fs(
    Y = Y_sum, 
    X = X_wo_1s, 
    k_constr = k_constr, 
    j_constr = 1, ## because it's been moved to the first column
    constraint_fn = constraint_fn,       
    constraint_grad_fn = constraint_grad_fn,
    maxit = maxit, 
    tol = tol,
    ls_max = ls_max, 
    ls_rho = ls_rho, 
    max_step = max_step,
    trackB = trackB,
    verbose = verbose
  )

  B <- out$B
  B <- B[, order(new_order)]
  
  # if we reassigned j_ref, adjust B so original reference is restored
  # if (!is.null(arch_j_ref)) {
  #   for (k in 1:nrow(B)) {
  #     B[k, ] <- B[k, ] - B[k, arch_j_ref]
  #   }
  # }
  
  z <- update_z(Y_orig, X, B)
  log_means <- X %*% B + matrix(z, ncol = 1) %*% matrix(1, ncol = J, nrow = 1)
  ll_new <- sum(Y_orig * log_means - exp(log_means))
  out$it_df[nrow(out$it_df), "lik"] <- ll_new
  
  ### now reorder B to match what was input 
  return(list(
    "B" = B,
    "k_constr" = k_constr,
    "j_constr" = j_constr,
    "niter" = out$iter,
    "loglik" = ll_new,
    "it_df" = out$it_df,
    "Bs" = out$history,
    "converged" = (out$iter < maxit)
  ))
  
}