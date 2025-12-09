#' fits model with B_kj constrained to equal g(B_k) for constraint fn g
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
#' @param constraint_tol tolerance for \eqn{\lvert B_kj - g(B_k)\rvert}
#' @param max_step maximum step size
#' @param c1 constant for armijo rule
#' @param maxit maximum iterations
#' @param inner_maxit max iterations per inner loop
#' @param verbose shout at you?
#' @param trackB track value of beta across iterations and return?
#' @param ignore_stop whether to ignore stopping criteria and run `maxit` iterations (could be helpful for diagnostic plots).
#' @param null_diagnostic_plots logical: should diagnostic plots be made for estimation under the null hypothesis? Default is \code{FALSE}.
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
#' equal to TRUE (otherwise it contains a NULL value).
#'
fit_null <- function(
  B,
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
  trackB = FALSE,
  ignore_stop = FALSE,
  null_diagnostic_plots = FALSE
) {
  J <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)

  it_df <- data.frame(
    it = 1:maxit,
    lik = NA,
    max_abs_B = NA,
    constraint_diff = NA
  )
  if (null_diagnostic_plots) {
    it_df$test_stat <- NA
  }

  #in vector format, which indexes to remove (bc B^j_constr = 0 id. constr.)
  indexes_to_remove <- (j_ref - 1) * p + 1:p

  #change id. constr. by subtracting approp. col of B from others
  for (k in 1:p) {
    B[k, ] <- B[k, ] - B[k, j_ref]
  }

  #update z
  z <- update_z(X = X, Y = Y, B = B)

  log_mean <- X %*% B + matrix(z, ncol = 1) %*% matrix(1, ncol = J, nrow = 1)

  ll <- sum(log_mean * Y - exp(log_mean))

  if (trackB) {
    Bs <- data.frame(k = rep(1:p, each = J), j = rep(1:J, p), B = NA)
  } else {
    Bs <- NULL
  }

  #get X_cup for later use
  if (is.null(X_cup)) {
    X_cup = X_cup_from_X(X, J)
  }

  #set iteration to zero
  iter <- 0

  converged <- FALSE

  #compute gap (i.e. g(B_k) - B_kj)
  gap <- constraint_fn[[k_constr]](B[k_constr, ]) - B[k_constr, j_constr]
  init_gap <- abs(gap)

  #set rho equal to initial value
  rho <- rho_init

  #initiate u
  u <- rho * gap

  if (trackB) {
    Bs$iter <- 0
    Bs$inner_iter <- 0
    Bs$rho <- rho
    Bs$u <- u

    Bs$gap <- gap
    #   Bs$score_stat <- 0
    lilBs <- Bs
    #
    for (k in 1:p) {
      Bs$B[1:J + J * (k - 1)] <- B[k, ]
    }
  }

  use_max_inner_it <- FALSE
  B_diff <- Inf

  criteria <- TRUE

  while (criteria) {
    iter <- iter + 1 #increment iteration
    inner_iter <- 0 #initiate internal iteration

    #evaluate augmented Lagrangian

    #  get current value of augmented lagrangian
    curr_lag_val <- -ll + u * gap + (rho / 2) * (gap^2)

    old_gap <- gap # old value of g(B_k) - B_kj
    old_B <- B #old value of B
    inner_diff <- Inf #initiate "observed" max abs val diff in B between inner iterations
    #at Inf so inner loop runs at least once
    #inner loop:
    while (
      ((inner_diff > inner_tol | use_max_inner_it) &
        inner_iter <= inner_maxit) |
        inner_iter == 0
    ) {
      inner_old_B <- B

      #perform block update to B as described in null estimation section of Clausen & Willis (2024)
      update <- macro_fisher_null(
        X = X,
        Y = Y,
        B = B,
        z = z,
        J = J,
        p = p,
        k_constr = k_constr,
        j_constr = j_constr,
        j_ref = j_ref,
        rho = rho,
        u = u,
        max_step = max_step,
        constraint_fn = constraint_fn,
        constraint_grad_fn = constraint_grad_fn,
        verbose = verbose,
        c1 = 1e-4
      )

      B <- B + update$update
      gap <- update$gap
      z <- update_z(Y, X, B)

      if (trackB) {
        temp_Bs <- lilBs
        for (k in 1:p) {
          temp_Bs$B[1:J + J * (k - 1)] <- B[k, ]
        }
        temp_Bs$iter <- iter
        temp_Bs$inner_iter <- inner_iter
        temp_Bs$rho <- rho
        temp_Bs$u <- u
        temp_Bs$gap <- gap
        Bs <- rbind(Bs, temp_Bs)
      }

      #did we move much?
      inner_diff <- max(abs(B - inner_old_B))

      #increment inner iteration counter
      inner_iter <- inner_iter + 1
    }
    #did we move much since last *outer* loop iteration?
    B_diff <- max(abs(B - old_B))

    if (verbose) {
      message(
        "Max absolute difference in B since last augmented Lagrangian outer step: ",
        signif(B_diff, 3)
      )
    }

    #update u and rho
    if (abs(gap) > constraint_tol) {
      u <- u + rho * gap
      if (abs(gap / old_gap) > kappa) {
        rho <- rho * tau
      }
    }

    if (verbose) {
      message("Estimate deviates from feasibility by ", signif(gap, 3))
      message("Parameter u set to ", signif(u, 3))
      message("Parameter rho set to ", signif(rho, 3))
    }

    if (abs(gap) < constraint_tol) {
      use_max_inner_it <- TRUE
    } else {
      use_max_inner_it <- FALSE
    }

    if (ignore_stop) {
      criteria <- (iter < maxit & !is.infinite(rho))
    } else {
      criteria <- (abs(gap) > constraint_tol | B_diff > B_tol) &
        iter < maxit &
        !is.infinite(rho)
    }

    log_means <- X %*% B + matrix(z, ncol = 1) %*% matrix(1, ncol = J, nrow = 1)

    ll <- sum(Y * log_means - exp(log_means))
    it_df$lik[iter] <- ll
    it_df$max_abs_B[iter] <- B_diff
    it_df$constraint_diff[iter] <- abs(gap)

    if (verbose) {
      message("Log likelihood: ", signif(ll, 3))
    }

    if (null_diagnostic_plots) {
      score_res <- try(get_score_stat(
        X = X,
        Y = Y,
        X_cup = X_cup,
        B = B,
        k_constr = k_constr,
        j_constr = j_constr,
        constraint_grad_fn = constraint_grad_fn,
        indexes_to_remove = (j_ref - 1) * p + 1:p,
        j_ref = j_ref,
        J = J,
        n = n,
        p = p
      ))
      if (!inherits(score_res, "try-error")) {
        it_df$test_stat[iter] <- score_res$score_stat
      }
    }

    if (abs(gap) <= constraint_tol & B_diff <= B_tol) {
      converged <- TRUE
    }
  }

  it_df <- it_df[1:iter, ]

  return(list(
    "B" = B,
    "k_constr" = k_constr,
    "j_constr" = j_constr,
    "niter" = iter,
    "gap" = gap,
    "u" = u,
    "rho" = rho,
    "Bs" = Bs,
    "it_df" = it_df,
    "converged" = converged
  ))
}
