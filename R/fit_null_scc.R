#' fits model with B_kj constrained to equal g(B_k) for constraint fn g, for a single category constraint
#' 
#' @param B description
#' @param Y Y (with augmentations)
#' @param X design matrix
#' @param X_cup design matrix for Y in long format. Defaults to NULL, in which case matrix is computed from X.
#' @param k_constr row index of B to constrain
#' @param j_constr col index of B to constrain
#' @param j_ref column index of convenience constraint
#' @param constraint_fn constraint function
#' @param constraint_grad_fn gradient of constraint fn
#' @param rho_init where to start quadratic penalty parameter
#' @param tau how much to increment rho by each iteration
#' @param kappa cutoff above which to increment rho. If distance to feasibility doesn't shrink by at least this factor in an iteration, increment rho by tau.
#' @param B_tol tolerance for convergence in \eqn{max_{k,j} \lvert B^t_{kj} - B^{(t - 1)}_{kj}\rvert}
#' @param inner_tol tolerance for inner loop
#' @param constraint_tol tolerance for \eqn{\lvert B_kj - g(B_k)\rvert }
#' @param max_step maximum step size
#' @param c1 constant for armijo rule
#' @param maxit maximum iterations
#' @param inner_maxit max iterations per inner loop
#' @param verbose shout at you?
#' @param trackB track value of beta across iterations and return?
#' 
#' @return A list containing elements `B`, `k_constr`, `j_constr`, `niter`
#' `gap`, `u`, `rho`, and `Bs`. `B` is a matrix containing parameter estimates 
#' under the null (obtained by maximum likelihood on augmented observations Y),
#' `k_constr`, and `j_constr` give row and column indexes of the parameter 
#' fixed to be equal to the constraint function \eqn{g()} under the null. `niter` is a 
#' scalar giving total number of outer iterations used to fit the null model, 
#' `gap` gives the final value of \eqn{g(B_{k constr}) - B_{k constr, j constr}}, 
#' `u` and `rho` are final values of augmented Lagrangian parameters, and 
#' `Bs` is a data frame containing values of B by iteration if `trackB` was set 
#' equal to TRUE (otherwise it contains a NULL value). - update based on new algorithm
#' 
fit_null_scc <- function(B,
                         Y, 
                         X, 
                         X_cup = NULL,
                         k_constr, 
                         j_constr, 
                         j_ref, 
                         constraint_fn, 
                         constraint_grad_fn, 
                         rho_init = 1, 
                         tau = 1.2, 
                         kappa = 0.8, 
                         B_tol = 1e-2, 
                         inner_tol = 0.01, 
                         constraint_tol = 1e-4, 
                         max_step = 5, 
                         c1 = 1e-4, 
                         maxit = 1000, 
                         inner_maxit = 25, 
                         verbose = FALSE, 
                         trackB = FALSE 
) {
  
  ### David to replace - for now redirecting to fit_null so that everything still runs! 
  ### AW Oct 31 2025: This functions looks redundant, but it's here because (I think)
  ### David said this could be sped up. So, this is really a placeholder for more complex code. 
  res <- fit_null(B = B,
                  Y = Y, 
                  X = X, 
                  X_cup = X_cup,
                  k_constr = k_constr, 
                  j_constr = j_constr, 
                  j_ref = j_ref, 
                  constraint_fn = constraint_fn, 
                  constraint_grad_fn = constraint_grad_fn, 
                  rho_init = rho_init, 
                  tau = tau, 
                  kappa = kappa, 
                  B_tol = B_tol, 
                  inner_tol = inner_tol, 
                  constraint_tol = constraint_tol, 
                  max_step = max_step, 
                  c1 = c1, 
                  maxit = maxit, 
                  inner_maxit = inner_maxit, 
                  verbose = verbose, 
                  trackB = trackB)
  B <- res$B
  iter <- res$niter
  gap <- res$gap
  u <- res$u
  rho <- res$rho
  Bs <- res$Bs
  ###
  
  return(list("B" = B,
              "k_constr" = k_constr,
              "j_constr" = j_constr,
              "niter" = iter,
              "gap" = gap,
              "u" = u,
              "rho" = rho,
              "Bs" = Bs))
}