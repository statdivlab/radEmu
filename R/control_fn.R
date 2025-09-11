#' Create list of control options (to pass to \code{emuFit()})
#'
#' @param control Current control list (optional), will augment it with missing arguments
#' @param max_step Maximum stepsize; update directions computed during estimation (under the alternative).
#' Will be rescaled if a step in any parameter exceeds this value. Defaults to 1.
#' @param ignore_stop whether to ignore stopping criteria and run `maxit` iterations (helpful for diagnostic plots to determine convergence).
#' @param use_fullmodel_info Used in estimation under the null hypothesis. Whether to use information matrix from
#' estimation under the alternative hypothesis to construct the robust score statistic (instead of information 
#' matrix from estimation under the null hypothesis). Defaults to `FALSE`.
#' @param use_fullmodel_cov Used in estimation under the null hypothesis. Whether to use covariance matrix from
#' estimation under the alternative hypothesis to construct the robust score statistic (instead of covariance 
#' matrix from estimation under the null hypothesis). Defaults to `FALSE`.
#' @param use_both_cov Used in estimation under the null hypothesis. Whether to do score test twice, once with 
#' covariance matrix under the alternative hypothesis and once with covariance matrix under the null hypothesis.
#' Defaults to `FALSE`. 
#' @param inner_maxit Used in estimation under the null hypothesis. Maximum number of iterations within each inner loop 
#' of estimation under null hypothesis algorithm. Default is `25`.
#' @param inner_tol Used in estimation under the null hypothesis. Convergence tolerance within each inner loop 
#' of estimation under null hypothesis algorithm. Default is `1`.
#' @param c1 Used in estimation under the null hypothesis. Parameter for Armijo line search. Default is `1e-4`.
#' @param trackB Used in estimation under the null hypothesis. When `TRUE` will track the value of `B` in each iteration
#' of optimization algorithm. Defaults to `FALSE`. 
#' @param return_nullB Used in estimation under the null hypothesis. When `TRUE` will return the final value of `B` 
#' under each null hypothesis tested. Defaults to `FALSE`. 
#' @param return_score_components Used in estimation under the null hypothesis. When `TRUE` will return the components
#' of the robust score test statistic for each null hypothesis tested. Defaults to `FALSE`. 
#' @param return_both_score_pvals Used in estimation under the null hypothesis, with `use_both_cov`. Defaults to `FALSE`. 
#' @param B_null_tol Used in estimation under the null hypothesis, for the augmented Lagrangian algorithm (by dafault
#' this algorithm is not used). numeric: convergence tolerance for null model fits for score testing (if max of absolute difference in
#' B across outer iterations is below this threshold, we declare convergence). Default is `0.001`.
#' @param rho_init Used in estimation under the null hypothesis, for the augmented Lagrangian algorithm (by default
#' this algorithm is not used). Value at which to initiate rho parameter in augmented Lagrangian
#' algorithm. Default is 1.  
#' @param tau Used in estimation under the null hypothesis, for the augmented Lagrangian algorithm (by default
#' this algorithm is not used). Value to scale `rho` by in each iteration of augmented Lagrangian
#' algorithm that does not move estimate toward zero sufficiently. Default is `2`.
#' @param kappa Used in estimation under the null hypothesis, for the augmented Lagrangian algorithm (by default
#' this algorithm is not used). Value between `0` and `1` that determines the cutoff on the ratio
#' of current distance from feasibility over distance in last iteration triggering
#' scaling of `rho`. If this ratio is above `kappa`, `rho` is scaled by `tau` to encourage
#' estimate to move toward feasibility.
#' @param constraint_tol Used in estimation under the null hypothesis, for the augmented Lagrangian algorithm (by default
#' this algorithm is not used). Constraint tolerance for fits under null hypotheses
#' (tested element of `B` must be equal to constraint function to within this tolerance for
#' a fit to be accepted as a solution to constrained optimization problem). Default is `1e-5`.  
#' @param ntries Used in estimation under the null hypothesis, for the augmented Lagrangian algorithm (by default
#' this algorithm is not used). The number of times to try optimization. Successive tries will change \code{tau} and
#' \code{inner_maxit} and retry.
#'
#' @return A list containing control options, to have more control over optimization algorithms used by `radEmu`. 
#' This can be passed into `emuFit()`.
#'
#' @export
control_fn <- function(
    control = list(), 
    max_step = 1,  
    ignore_stop = FALSE,
    use_fullmodel_info = FALSE,
    use_fullmodel_cov = FALSE,
    use_both_cov = FALSE,
    inner_maxit = 25, 
    inner_tol = 1, 
    c1 = 1e-4,
    trackB = FALSE,
    return_nullB = FALSE,
    return_score_components = FALSE,
    return_both_score_pvals = FALSE,
    B_null_tol = 0.001,
    rho_init = 1,
    tau = 2,
    kappa = 0.8,
    constraint_tol = 1e-5,
    ntries = 4
    ) {
  
  default <- list(max_step = max_step,
                  ignore_stop = ignore_stop, 
                  use_fullmodel_info = use_fullmodel_info,
                  use_fullmodel_cov = use_fullmodel_cov,
                  use_both_cov = use_both_cov,
                  inner_maxit = inner_maxit, 
                  inner_tol = inner_tol, 
                  c1 = c1,
                  trackB = trackB,
                  return_nullB = return_nullB,
                  return_score_components = return_score_components,
                  return_both_score_pvals = return_both_score_pvals,
                  B_null_tol = B_null_tol,
                  rho_init = rho_init,
                  tau = tau,
                  kappa = kappa,
                  constraint_tol = constraint_tol,
                  ntries = ntries)
  
  if (length(control) > 0) {
    unknown <- setdiff(names(control), names(default))
    if (length(unknown) > 0) {
      warning("Unknown control parameters: ", paste(unknown, collapse = ", "))
    }
    default[names(control)] <- control
  }
  
  return(default)
}