#' Run robust score test

#' @param B value of coefficient matrix (p x J) returned by full model fit or value of coefficient
#' matrix to start null estimation at given as input to emuFit
#' @param Y an n x J matrix or dataframe of *augmented* nonnegative observations (i.e.,
#' observations Y plus augmentations from last iteration of maximum penalized likelihood estimation
#' for full model)
#' @param X an n x p matrix of covariates
#' @param X_cup the design matrix for long format Y in long format B (nJ x pJ)
#' @param k_constr row index of element of B to be tested for equality to row identifiability constraint
#' @param j_constr column index of element of B to be tested for equality to row identifiability constraint
#' @param constraint_fn constraint function g to be imposed on each row of B; the null that
#' is tested is $B_{k_constr,j_constr} = g(B_k_constr)$
#' @param constraint_grad_fn derivative of constraint_fn with respect to its
#' arguments (i.e., elements of a row of B)
#' @param j_ref column dropped from B in convenience parametrization used to compute
#' components of score statistic and fit null; test is invariant to choice of j_ref,
#' but optimization tends to be faster and more stable if j_ref is a column of Y with
#' many positive entries
#' @param c1 numeric: parameter for Armijo line search. Default is 1e-4.
#' @param maxit numeric: maximum number of outer iterations to perform (default is 1000)
#' @param inner_maxit numeric: maximum number of inner iterations to perform
#' per outer iteration (default is 25)
#' is not provided, all elements of B save the intercept row will be tested.
#' @param constraint_fn function g defining a constraint on rows of B; g(B_k) = 0
#' for rows k = 1, ..., p of B. Default function is a smoothed median (minimizer of
#' pseudohuber loss).
#' @param constraint_grad_fn function returning gradient of constraint function (as
#' a function of a row of B)
#' @param rho_init numeric: value at which to initiate rho parameter in augmented Lagrangian
#' algorithm. Default is 1.
#' @param tau numeric: value to scale rho by in each iteration of augmented Lagrangian
#' algorithm that does not move estimate toward zero sufficiently. Default is 2.
#' @param kappa numeric: value between 0 and 1 that determines the cutoff on the ratio
#' of current distance from feasibility over distance in last iteration triggering
#' scaling of rho. If this ratio is above kappa, rho is scaled by tau to encourage
#' estimate to move toward feasibility.
#' @param B_tol numeric: convergence tolerance for null model fits for score testing (if max of absolute difference in
#' B across outer iterations is below this threshold, we declare convergence).
#' Default is 0.001.
#' @param inner_tol numeric: convergence tolerance for inner loop of null fitting algorithm (if max of absolute difference in
#' B across inner iterations is below this threshold, we declare convergence).
#' Default is 0.01.
#' @param constraint_tol numeric: constraint tolerance for fits under null hypotheses
#' (tested element of B must be equal to constraint function to within this tolerance for
#' a fit to be accepted as a solution to constrained optimization problem). Default is 1e-5.
#' @param j_ref column index of convenience constraint
#' @param c1 numeric: parameter for Armijo line search. Default is 1e-4.
#' @param maxit maximum number of outer iterations of augmented lagrangian algorithm to perform before
#' exiting optimization. Default is 1000.
#' @param inner_maxit maximum number of coordinate descent passes through columns of B to make within each
#' outer iteration of augmented lagrangian algorithm before exiting inner loop
#' @param ntries numeric: total number of times to attempt optimization under null
#' if optimization fails (optimization parameters will be tweaked in subsequent
#' fits to attempt to avoid failure). Default is 4.
#' @param verbose provide updates as model is being fitted? Defaults to TRUE.
#' @param trackB store and return values of B at each iteration of optimization
#' algorithm? Useful for debugging. Default is FALSE.
#' @param I_inv Optional: matrix containing inverted information matrix computed
#' under full model. Default is NULL, in which case information is recomputed under
#' null, which we recommend.
#' @param Dy Optional: matrix containing empirical score covariance computed
#' under full model. Default is NULL, in which case this quantity is recomputed
#' under null, which we recommend.
#' @param return_both_score_pvals logical: should score p-values be returned using both
#' information matrix computed from full model fit and from null model fits? Default is
#' FALSE. This parameter is used for simulations - in any applied analysis, type of
#' p-value to be used should be chosen before conducting tests.
#' @param cluster a numeric vector giving cluster membership for each row of Y to
#' be used in computing GEE test statistics. Default is NULL, in which case rows of
#' Y are treated as independent.
#' @param null_fit_constraint the type of constraint, which informs which algorithm
#' will be used to fit the model under the null hypothesis. NULL by default, in which
#' case the standard fitting algorithm will be used. If included, this argument must be either
#' `scc` for single category constraint, `symmetric` for symmetric functions such
#' as mean or pseudo-Huber median, or `other`.
#'
#' @return A list containing elements `score_stat`, `pval`, `log_pval`,'niter`,
#' `convergence`, `gap`, `u`, `rho`, `tau`, `inner_maxit`, `null_B`, and `Bs`. `score_stat` gives the
#' value of the robust score statistic for $H_0: B_{k_constr,j_constr} = g(B_{k_constr})$.
#' `pval` and `log_pval` are the p-value (on natural and log scales) corresponding to
#' the score statistic (log_pval may be useful when the p-value is very close to zero).
#' `gap` is the final value of $g(B_{k_constr}) - B_{k_constr, j_constr}$ obtained in
#' optimization under the null. `u` and `rho` are final values of augmented
#' Lagrangian parameters returned by null fitting algorithm. `tau` is the final value of `tau` that
#' is used to update the `rho` values and `inner_maxit` is the final maximum number of iterations for
#' the inner optimization loop in optimization under the null, in which B and z parameter values are
#' maximized for specific `u` and `rho` parameters. `null_B` is the value of
#' B returned but the null fitting algorithm. `Bs` is by default `NULL`; if `trackB = TRUE`,
#' `Bs` is a data frame containing values of B by outcome category, covariate, and
#' iteration.
#'
#'
#'
#' @importFrom stats cov median model.matrix optim pchisq qnorm weighted.mean
#' @import Matrix
#' @import MASS
#'
#' @export
score_test <- function(
  B, #B (MPLE)
  Y, #Y (with augmentations)
  X, #design matrix
  X_cup,
  k_constr, #row index of B to constrain
  j_constr,
  constraint_fn, #constraint function
  constraint_grad_fn, #gradient of constraint fn
  rho_init = 1,
  tau = 2,
  kappa = 0.8,
  B_tol = 1e-3,
  inner_tol = 0.01,
  constraint_tol = 1e-3,
  j_ref,
  c1 = 1e-4,
  maxit = 1000,
  inner_maxit = 25,
  ntries = 4,
  verbose = FALSE,
  trackB = FALSE,
  I_inv = NULL,
  Dy = NULL,
  return_both_score_pvals = FALSE,
  cluster = NULL,
  null_fit_constraint = NULL
) {
  # get hyperparameters
  n <- nrow(Y)
  J <- ncol(Y)
  p <- ncol(X)

  # identify if constraint function is symmetric or a single category constraint
  if (is.null(null_fit_constraint)) {
    constraint_type <- "other"
  } else {
    constraint_type <- null_fit_constraint
  }
  if (!(constraint_type %in% c("other", "symmetric", "scc"))) {
    stop(
      "The argument `null_fit_constraint` must either be omitted, or must be either
         `scc` for single category constraint, `symmetric` for symmetric functions such
         as mean or pseudo-Huber median, or `other`."
    )
  }

  # start fitting
  tries_so_far <- 0
  accept_try <- FALSE
  good_enough_fit <- FALSE
  while (!accept_try) {
    #fit under null
    # note that the below could be more concise by dynamically selecting function based on constraint
    # type without the if/else, however that would break if these functions needed different arguments,
    # which seemed possible. If they end up having the same arguments, Sarah will fix this
    if (constraint_type == "other") {
      constrained_fit <- try(fit_null(
        B = B, #B (MPLE)
        Y = Y, #Y (with augmentations)
        X = X, #design matrix
        X_cup = X_cup,
        k_constr = k_constr, #row index of B to constrain
        j_constr = j_constr, #col index of B to constrain
        constraint_fn = constraint_fn, #constraint function
        constraint_grad_fn = constraint_grad_fn, #gradient of constraint fn
        rho_init = rho_init,
        tau = tau,
        kappa = kappa,
        B_tol = B_tol,
        inner_tol = inner_tol,
        constraint_tol = constraint_tol,
        j_ref = j_ref,
        c1 = c1,
        maxit = maxit,
        inner_maxit = inner_maxit,
        verbose = verbose,
        trackB = trackB
      ))
    } else if (constraint_type == "scc") {
      constrained_fit <- try(fit_null_scc(
        B = B, #B (MPLE)
        Y = Y, #Y (with augmentations)
        X = X, #design matrix
        X_cup = X_cup,
        k_constr = k_constr, #row index of B to constrain
        j_constr = j_constr, #col index of B to constrain
        constraint_fn = constraint_fn, #constraint function
        constraint_grad_fn = constraint_grad_fn, #gradient of constraint fn
        rho_init = rho_init,
        tau = tau,
        kappa = kappa,
        B_tol = B_tol,
        inner_tol = inner_tol,
        constraint_tol = constraint_tol,
        j_ref = j_ref,
        c1 = c1,
        maxit = maxit,
        inner_maxit = inner_maxit,
        verbose = verbose,
        trackB = trackB
      ))
    } else {
      constrained_fit <- try(fit_null_symmetric(
        B = B, #B (MPLE)
        Y = Y, #Y (with augmentations)
        X = X, #design matrix
        X_cup = X_cup,
        k_constr = k_constr, #row index of B to constrain
        j_constr = j_constr, #col index of B to constrain
        constraint_fn = constraint_fn, #constraint function
        constraint_grad_fn = constraint_grad_fn, #gradient of constraint fn
        B_tol = B_tol,
        j_ref = j_ref,
        c1 = c1,
        maxit = maxit,
        inner_maxit = inner_maxit,
        verbose = verbose,
        trackB = trackB
      ))
    }

    if (inherits(constrained_fit, "try-error")) {
      accept_try <- FALSE
    } else {
      if (
        (abs(constrained_fit$gap) <= constraint_tol) &
          (constrained_fit$niter < maxit) &
          (!is.infinite(constrained_fit$rho))
      ) {
        accept_try <- TRUE
        good_enough_fit <- TRUE
      } else {
        tau <- tau^(3 / 4)
        inner_maxit <- 2 * inner_maxit
        message(
          "Constrained optimization failed to converge within iteration limit;
retrying with smaller penalty scaling parameter tau and larger inner_maxit."
        )
      }
    }
    tries_so_far <- tries_so_far + 1
    if (tries_so_far == ntries) {
      accept_try <- TRUE
    }
  }

  if (!good_enough_fit & constraint_type == "symmetric") {
    if (inherits(constrained_fit, "try-error")) {
      message("Optimization via Fisher scoring failed, attempting augmented Lagrangian approach.")
      constrained_fit <- try(fit_null(
        B = B, #B (MPLE)
        Y = Y, #Y (with augmentations)
        X = X, #design matrix
        X_cup = X_cup,
        k_constr = k_constr, #row index of B to constrain
        j_constr = j_constr, #col index of B to constrain
        constraint_fn = constraint_fn, #constraint function
        constraint_grad_fn = constraint_grad_fn, #gradient of constraint fn
        rho_init = rho_init,
        tau = tau,
        kappa = kappa,
        B_tol = B_tol,
        inner_tol = inner_tol,
        constraint_tol = constraint_tol,
        j_ref = j_ref,
        c1 = c1,
        maxit = maxit,
        inner_maxit = inner_maxit,
        verbose = verbose,
        trackB = trackB
      ))
      
      if (!inherits(constrained_fit, "try-error")) {
        if (
          (abs(constrained_fit$gap) <= constraint_tol) &
          (constrained_fit$niter < maxit) &
          (!is.infinite(constrained_fit$rho))
        ) {
          accept_try <- TRUE
          good_enough_fit <- TRUE
        }
      }
    }
  }
  
  if (!good_enough_fit) {
    
    warning(
      "Optimization for null fit with k = ",
      k_constr,
      " and j = ",
      j_constr,
      " failed to converge across ",
      ntries,
      ifelse(ntries > 1, " attempts.", " attempt.")
    )
    
  }
  B <- constrained_fit$B
  z <- update_z(Y, X, B)
  p <- ncol(X)
  J <- ncol(Y)
  n <- nrow(Y)

  # message("Make X_cup an argument of score_test and calculate it once in emuFit")
  # X_cup = X_cup_from_X(X,J)
  #indexes in long format corresponding to the j_constr-th col of B
  #get score stat
  indexes_to_remove <- (j_ref - 1) * p + 1:p
  score_res <- try(
    get_score_stat(
      Y = Y,
      X_cup = X_cup,
      X = X,
      B = B,
      k_constr = k_constr,
      j_constr = j_constr,
      j_ref = j_ref,
      constraint_grad_fn = constraint_grad_fn,
      indexes_to_remove = indexes_to_remove,
      J = J,
      n = n,
      p = p,
      I_inv = I_inv,
      Dy = Dy,
      cluster = cluster
    )
  )
  if (inherits(score_res, "try-error")) {
    score_stat <- score_res
  } else {
    score_stat <- score_res$score_stat
  }

  if (!return_both_score_pvals) {
    #typically we want only one score p-value
    #(using only one version of information matrix)
    if (inherits(score_stat, "try-error")) {
      warning(
        "score statistic for test of k = ",
        k_constr,
        " and j = ",
        j_constr,
        " cannot be computed, likely because the information matrix is computationally singular."
      )
      score_stat <- NA
    }
    return(list(
      "score_stat" = score_stat,
      "score_pieces" = score_res,
      "pval" = pchisq(score_stat, 1, lower.tail = FALSE),
      "log_pval" = pchisq(score_stat, 1, lower.tail = FALSE, log.p = TRUE),
      "niter" = constrained_fit$niter,
      "convergence" = ifelse(
        constrained_fit$niter >= maxit,
        'iteration limit reached',
        ifelse(
          is.infinite(constrained_fit$rho),
          'problem became ill-conditioned',
          'converged'
        )
      ),
      # "proj_score" = constrained_fit$proj_score,
      "gap" = constrained_fit$gap,
      "u" = constrained_fit$u,
      "rho" = constrained_fit$rho,
      "tau" = tau,
      "inner_maxit" = inner_maxit,
      "null_B" = constrained_fit$B,
      # "score_stats" = constrained_fit$score_stats,
      "Bs" = constrained_fit$Bs,
      "niter" = constrained_fit$niter
    ))
  } else {
    #for simulations -- if we want to return both the score p-value using
    #information from full model fit and from null model
    score_res_with_null_info <-
      get_score_stat(
        Y = Y,
        X_cup = X_cup,
        X = X,
        B = B,
        k_constr = k_constr,
        j_constr = j_constr,
        j_ref = j_ref,
        constraint_grad_fn = constraint_grad_fn,
        indexes_to_remove = indexes_to_remove,
        J = J,
        n = n,
        p = p,
        I_inv = NULL,
        Dy = Dy
      )
    if (inherits(score_res_with_null_info, "try-error")) {
      score_stat_with_null_info <- score_res_with_null_info
    } else {
      score_stat_with_null_info <- score_res_with_null_info$score_stat
    }
    score_stat_with_null_info <- score_stat_with_null_info
    if (inherits(score_stat_with_null_info, "try-error")) {
      warning(
        "one of the score statistics for test of k = ",
        k_constr,
        " and j = ",
        j_constr,
        " cannot be computed, likely because the information matrix is computationally singular."
      )
      score_stat_with_null_info <- NA
    }

    return(list(
      "score_stat" = score_stat,
      "score_pieces" = score_res,
      "pval" = pchisq(score_stat, 1, lower.tail = FALSE),
      "log_pval" = pchisq(score_stat, 1, lower.tail = FALSE, log.p = TRUE),
      "score_stat_null_info" = score_stat_with_null_info,
      "score_pieces_null_info" = score_res_with_null_info,
      "pval_null_info" = pchisq(
        score_stat_with_null_info,
        1,
        lower.tail = FALSE
      ),
      "log_pval_null_info" = pchisq(
        score_stat_with_null_info,
        1,
        lower.tail = FALSE,
        log.p = TRUE
      ),
      "niter" = constrained_fit$niter,
      "convergence" = ifelse(
        constrained_fit$niter >= maxit,
        'iteration limit reached',
        'converged'
      ),
      # "proj_score" = constrained_fit$proj_score,
      "gap" = constrained_fit$gap,
      "u" = constrained_fit$u,
      "rho" = constrained_fit$rho,
      "tau" = tau,
      "inner_maxit" = inner_maxit,
      "null_B" = constrained_fit$B,
      # "score_stats" = constrained_fit$score_stats,
      "Bs" = constrained_fit$Bs,
      "niter" = constrained_fit$niter
    ))
  }
}
