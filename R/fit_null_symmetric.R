#' fits model with B_kj constrained to equal g(B_k) for constraint fn g, for a symmetric constraint
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
#' @param B_tol tolerance for convergence in $max_{k,j} |B^t_{kj} - B^{(t - 1)}_{kj}|$
#' @param inner_tol tolerance for inner loop
#' @param c1 constant for armijo rule
#' @param maxit maximum iterations
#' @param inner_maxit max iterations per inner loop
#' @param verbose shout at you?
#' @param trackB track value of beta across iterations and return?
#' @param use_optim whether to use `optim` instead of fisher scoring. Default is FALSE.
#'
#' @return A list containing elements `B`, `k_constr`, `j_constr`, `niter`
#' `gap`, `u`, `rho`, and `Bs`. `B` is a matrix containing parameter estimates
#' under the null (obtained by maximum likelihood on augmented observations Y),
#' `k_constr`, and `j_constr` give row and column indexes of the parameter
#' fixed to be equal to the constraint function $g()$ under the null. `niter` is a
#' scalar giving total number of outer iterations used to fit the null model,
#' `gap` gives the final value of $g(B_{k constr}) - B_{k constr, j constr}$,
#' `u` and `rho` are final values of augmented Lagrangian parameters, and
#' `Bs` is a data frame containing values of B by iteration if `trackB` was set
#' equal to TRUE (otherwise it contains a NULL value). - update based on new algorithm
#'
fit_null_symmetric <- function(
  B,
  Y,
  X,
  X_cup = NULL,
  k_constr,
  j_constr,
  j_ref,
  constraint_fn,
  constraint_grad_fn,
  B_tol = 1e-2,
  c1 = 1e-2,
  maxit = 1000,
  inner_maxit = 25,
  inner_tol = 0.01,
  verbose = FALSE,
  trackB = FALSE,
  use_optim = FALSE
) {
  if (is.list(constraint_fn)) {
    constraint_fn <- constraint_fn[[k_constr]]
    constraint_grad_fn <- constraint_grad_fn[[k_constr]]
  }
  #useful constants
  p <- ncol(X)
  n <- nrow(Y)
  J <- ncol(Y)

  #reparametrize if convenience constraint taxon is also taxon containing null-constrained param
  if (j_ref == j_constr) {
    arch_j_ref <- j_ref
    new_j_ref <- which.max(colSums(Y > 1)[-j_ref])
    new_j_ref <- new_j_ref + as.numeric(new_j_ref >= j_ref)
    j_ref <- new_j_ref
  } else {
    arch_j_ref <- NULL
  }

  #impose convenience constraint
  for (k in 1:p) {
    B[k, ] <- B[k, ] - B[k, j_ref]
  }

  #force B to satisfy null constraint
  #(this is the bit where symmetry of g() is important)
  B[k_constr, j_constr] <- constraint_fn(B[k_constr, -j_constr])

  #update z
  z <- update_z(Y, X, B)

  #create df for tracking B through optimization if trackB = TRUE
  if (trackB) {
    df <- data.frame(
      iter = numeric(0),
      iter_j = numeric(0),
      k = numeric(0),
      j = numeric(0),
      value = numeric(0),
      ll = numeric(0),
      gr_norm = numeric(0)
    )
  }

  #break taxa, j_ref and j_constr excluded, into two groups
  #according to whether they are above or below B[k_constr,j_constr]
  #and order these groups so each set begins with j s.t. B[k_constr,j] is
  #closest to B[k_constr,j_constr]

  above_g <- setdiff(
    which(B[k_constr, ] > constraint_fn(B[k_constr, ])),
    c(j_ref, j_constr)
  )
  above_g <- above_g[order(B[k_constr, above_g])]
  below_g <- setdiff(
    which(B[k_constr, ] < constraint_fn(B[k_constr, ])),
    c(j_ref, j_constr)
  )
  below_g <- below_g[order(-B[k_constr, below_g])]

  #count total pairings between taxa in above_g and below_g, and count remaining taxa
  npairs <- min(length(above_g), length(below_g))
  nsingles <- max(length(above_g), length(below_g)) - npairs

  #now create pairs sandwiching B[k_constr,j_constr] and put remaining taxa into 'singles'
  if (npairs > 0) {
    pairs <- cbind(above_g[1:npairs], below_g[1:npairs])
    singles <- c(above_g[-(1:npairs)], below_g[-(1:npairs)])
  } else {
    pairs <- integer(0)
    singles <- c(above_g, below_g)
  }

  #compute fitted log means
  log_means <- X %*% B
  for (i in 1:n) {
    log_means[i, ] <- log_means[i, ] + z[i]
  }

  #compute ll
  ll <- sum(Y * log_means - exp(log_means))

  #initialize 'keep_going' indicator for continuing optimization as well as 'iter'
  #to track optimization iterations
  keep_going <- TRUE
  iter <- 1

  # temporarily set use_optim to `TRUE` always
  # use_optim <- TRUE

  #iterate until something happens
  while (keep_going) {
    prev_B <- B
    prev_ll <- ll

    #loop through paired above_g and below_g taxa, and then through remaining singletons
    for (jset in 1:(npairs + nsingles)) {
      #choose groups of taxa to optimize over (in addition to j_constr, which is always included)
      if (jset > npairs) {
        js <- singles[jset - npairs]
      } else {
        js <- pairs[jset, ]
      }

      njs <- length(js)

      #computed update for js and j_constr
      if (use_optim) {
        update <- optim(
          rep(0, p * njs + p - 1),
          fn = function(x) {
            -null_repar_ll(
              x = x,
              js = js,
              B = B,
              z = z,
              p = p,
              Y = Y,
              X = X,
              j_constr = j_constr,
              k_constr = k_constr,
              constraint_fn = constraint_fn
            )
          },
          # gr = function(x) -null_repar_ll_gr(x = x,
          #                                    js = js,
          #                                    B = B,
          #                                    z = z,
          #                                    p = p,
          #                                    Y = Y,
          #                                    X = X,
          #                                    j_constr = j_constr,
          #                                    k_constr = k_constr,
          #                                    constraint_fn = constraint_fn,
          #                                    constraint_grad_fn = constraint_grad_fn),
          method = "L-BFGS-B",
          control = list(maxit = maxit),
          lower = -1,
          upper = 1
        )

        #unfortunately optim sometimes returns 'optimal' parameter values that
        #decrease the likelihood relative to starting parameters, so we have to
        #check that this has not occurred
        #(why, optim?)
        if (
          -update$value >
            null_repar_ll(
              x = rep(0, p * njs + p - 1),
              js = js,
              B = B,
              z = z,
              p = p,
              Y = Y,
              X = X,
              j_constr = j_constr,
              k_constr = k_constr,
              constraint_fn = constraint_fn
            )
        ) {
          #update elements of B corresponding to j in js
          for (jind in 1:njs) {
            B[, js[jind]] <- B[, js[jind]] + update$par[1:p + (jind - 1) * p]
          }
          #update unconstrained elements of B[,j_constr]
          B[-k_constr, j_constr] <- B[-k_constr, j_constr] +
            update$par[njs * p + 1:(p - 1)]
          #and update constrained element
          B[k_constr, j_constr] <- constraint_fn(B[k_constr, -j_constr])
        } else {
          if (verbose) {
            message(
              "Update for j = ",
              pairs[npairs, 1],
              " and ",
              pairs[npairs, 2],
              " did not increase log likelihood. Update skipped."
            )
          }
        }
      } else {
        # Fisher scoring (by Amy & copilot)

        #------------------------------------------------------
        # 1.  Build problem-specific closures once
        # -----------------------------------------------------
        make_fs_problem <- function(
          js,
          B,
          z,
          p,
          Y,
          X,
          j_constr,
          k_constr,
          constraint_fn,
          constraint_grad_fn
        ) {
          fn <- function(theta) {
            -null_repar_ll(
              x = theta,
              js = js,
              B = B,
              z = z,
              p = p,
              Y = Y,
              X = X,
              j_constr = j_constr,
              k_constr = k_constr,
              constraint_fn = constraint_fn
            )
          }

          grad_fn <- function(theta) {
            ## need to negate gradient later; can't do it here as fn returns a list
            null_repar_ll_gr(
              x = theta,
              js = js,
              B = B,
              z = z,
              p = p,
              Y = Y,
              X = X,
              j_constr = j_constr,
              k_constr = k_constr,
              constraint_fn = constraint_fn,
              constraint_grad_fn = constraint_grad_fn,
              return_hess = TRUE, # we want info + grad
              return_info_inv = FALSE # we don't want to invert info in grad fn
            )
          }

          list(fn = fn, grad_fn = grad_fn)
        }

        # -------------------------------------------------------
        # 2.  Generic Fisher-scoring optimiser
        # -------------------------------------------------------
        fisher_scoring <- function(
          theta0,
          fn,
          grad_fn,
          max_iter,
          c1, # Armijo constant
          backtrack_max,
          tol,
          verbose
        ) {
          theta <- theta0

          for (iter in seq_len(max_iter)) {
            # -- gradient and expected Fisher information
            g_obj <- grad_fn(theta)
            grad <- -as.numeric(g_obj$gr) ### negate here
            info <- g_obj$info

            # -- ensure info is invertible (lambda·diag trick)
            lambda <- 0
            solve_ok <- FALSE
            repeat {
              if (lambda > 0) {
                info_reg <- info +
                  lambda *
                    Matrix::Diagonal(
                      n = nrow(info),
                      x = pmax(abs(Matrix::diag(info)), 1)
                    )
              } else {
                info_reg <- info
              }

              step_dir <- try(Matrix::solve(info_reg, grad), silent = TRUE)
              if (!inherits(step_dir, "try-error")) {
                step_dir <- -as.numeric(step_dir) # Newton-style step (downhill)
                solve_ok <- TRUE
                break
              }
              lambda <- if (lambda == 0) 1e-4 else 10 * lambda
              if (lambda > 1e6) {
                stop("Unable to regularise Fisher information for inversion")
              }
            }

            # -- back-tracking Armijo line search
            step <- 1
            f_curr <- fn(theta)
            slope <- sum(grad * step_dir)

            for (bt in seq_len(backtrack_max)) {
              theta_new <- theta + step * step_dir
              f_new <- fn(theta_new)
              if (f_new <= f_curr + c1 * step * slope) {
                break
              }
              step <- step / 2
            }

            # -- update
            theta <- theta_new
            if (verbose) {
              cat(sprintf(
                "Iter %2d  f = %-12.6g  abs theta change = %-10.3g  lambda = %g\n",
                iter,
                f_new,
                max(abs(step * step_dir)),
                lambda
              ))
            }

            # -- convergence?
            if (max(abs(step * step_dir)) < tol) {
              return(list(
                par = theta,
                value = f_new,
                converged = TRUE,
                iter = iter
              ))
            }
          }

          list(
            par = theta,
            value = fn(theta),
            converged = FALSE,
            iter = max_iter
          )
        }

        # --------------------------------------------------------
        # 3.  Do it
        # --------------------------------------------------------
        prob <- make_fs_problem(
          js,
          B,
          z,
          p,
          Y,
          X,
          j_constr,
          k_constr,
          constraint_fn,
          constraint_grad_fn
        )

        update <- fisher_scoring(
          rep(0, p * length(js) + p - 1),
          fn = prob$fn,
          grad_fn = prob$grad_fn,
          #max_iter = maxit,
          max_iter = inner_maxit,
          c1 = c1, # Armijo constant
          backtrack_max = inner_maxit,
          #tol = B_tol,
          tol = inner_tol, 
          verbose = verbose
        )

        for (jind in 1:njs) {
          B[, js[jind]] <- B[, js[jind]] + update$par[1:p + (jind - 1) * p]
        }
        #update unconstrained elements of B[,j_constr]
        B[-k_constr, j_constr] <- B[-k_constr, j_constr] +
          update$par[njs * p + 1:(p - 1)]
        #and update constrained element
        B[k_constr, j_constr] <- constraint_fn(B[k_constr, -j_constr])
      }

      if (trackB) {
        gr <- get_constrained_gr(
          Y = Y,
          X = X,
          B = B,
          z = z,
          js = setdiff(1:J, c(j_constr, j_ref)),
          j_constr = j_constr,
          k_constr = k_constr,
          j_ref = j_ref,
          constraint_grad_fn = constraint_grad_fn,
          gr_only = TRUE
        )

        df <- rbind(
          df,
          data.frame(
            iter = iter,
            k = rep(1:p, each = length(js) + 1),
            j = rep(c(js, j_constr), p),
            value = do.call(c, lapply(1:p, function(k) B[k, c(js, j_constr)])),
            ll = -ll,
            gr_norm = sum(gr^2)
          )
        )
      }

      #update z
      z <- update_z(Y = Y, X = X, B = B)
    }

    #compute ll again
    log_means <- X %*% B
    for (i in 1:nrow(Y)) {
      log_means[i, ] <- log_means[i, ] + z[i]
    }

    ll <- sum(Y * log_means - exp(log_means))

    #compute gradient of ll

    if (verbose) {
      #compute gradient if we didn't already
      if (!trackB) {
        gr <- get_constrained_gr(
          Y = Y,
          X = X,
          B = B,
          z = z,
          js = setdiff(1:J, c(j_constr, j_ref)),
          j_constr = j_constr,
          k_constr = k_constr,
          j_ref = j_ref,
          constraint_grad_fn = constraint_grad_fn,
          gr_only = TRUE
        )
      }
      #tell you all about the ll, gradient, and how much B has changed this loop
      message("ll = ", round(ll, 1))
      message("ll increased by ", round(ll - prev_ll, 1))
      message("gr_norm = ", round(sum(gr^2), 1))
      message("max abs diff in B = ", signif(max(abs(B - prev_B)), 2))
    }

    #recompute above_g and below_g
    above_g <- setdiff(
      which(B[k_constr, ] > constraint_fn(B[k_constr, ])),
      c(j_ref, j_constr)
    )
    above_g <- above_g[order(B[k_constr, above_g])]
    below_g <- setdiff(
      which(B[k_constr, ] < constraint_fn(B[k_constr, ])),
      c(j_ref, j_constr)
    )
    below_g <- below_g[order(-B[k_constr, below_g])]

    #and pairs of taxa and singleton taxa to optimize over iteratively
    npairs <- min(length(above_g), length(below_g))
    nsingles <- max(length(above_g), length(below_g)) - npairs

    if (npairs > 0) {
      pairs <- cbind(above_g[1:npairs], below_g[1:npairs])
      singles <- c(above_g[-(1:npairs)], below_g[-(1:npairs)])
    } else {
      pairs <- integer(0)
      singles <- c(above_g, below_g)
    }

    #increment iteration counter
    iter <- iter + 1

    #determine whether we should stop iterating
    if (iter > maxit) {
      keep_going <- FALSE
      converged <- FALSE
    }

    if (max(abs(B - prev_B)) < B_tol) {
      keep_going <- FALSE
      converged <- TRUE
    }
  }

  #if we changed j_ref, change it back
  if (!is.null(arch_j_ref)) {
    for (k in 1:p) {
      B[k, ] <- B[k, ] - B[k, arch_j_ref]
    }
  }

  #create something to return if trackB is FALSE
  if (!trackB) {
    df = NULL
  }
  #return a list containing B 'n' fixin's
  return(list(
    B = B,
    k_constr = k_constr,
    j_constr = j_constr,
    niter = iter,
    converged = converged,
    Bs = df,
    gap = 0, #gap is zero by parametrization
    "u" = NA,
    "rho" = NA
  ))
}
