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
#' @param B_tol tolerance for convergence in \eqn{max_{k,j} \lvert B^t_{kj} - B^{(t - 1)}_{kj}\rvert}
#' @param inner_tol tolerance for inner loop
#' @param c1 constant for armijo rule
#' @param maxit maximum iterations
#' @param inner_maxit max iterations per inner loop
#' @param verbose shout at you?
#' @param trackB track value of beta across iterations and return?
#' @param use_optim whether to use `optim` instead of fisher scoring. Default is FALSE.
#' @param ignore_stop whether to ignore stopping criteria and run `maxit` iterations (could be helpful for diagnostic plots).
#' @param tol_lik tolerance for relative changes in likelihood for stopping criteria. Default is `1e-5`.
#' @param tol_test_stat tolerance for relative changes in test statistic for stopping criteria. Default is `0.01`.
#' @param null_window window to use for stopping criteria (this many iterations where stopping criteria is met). Default is `5`.
#' @param max_step Default is `1`.
#' @param reference_set Optional reference set to act as the subset. Default is `NULL` in which case the function will look into 
#' `constraint_fn` for the reference set.
#'
#' @return A list containing elements `B`, `k_constr`, `j_constr`, `niter`
#' and `Bs`. `B` is a matrix containing parameter estimates
#' under the null (obtained by maximum likelihood on augmented observations Y),
#' `k_constr`, and `j_constr` give row and column indexes of the parameter
#' fixed to be equal to the constraint function \eqn{g()} under the null. `niter` is a
#' scalar giving total number of outer iterations used to fit the null model, and
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
  inner_tol = 0.01,
  c1 = 1e-2,
  maxit = 1000,
  inner_maxit = 25,
  verbose = FALSE,
  trackB = FALSE,
  use_optim = FALSE,
  ignore_stop = FALSE,
  tol_lik = 1e-5,
  tol_test_stat = 0.01,
  null_window = 5,
  max_step = 1,
  reference_set = NULL
) {
  
  #useful constants
  p <- ncol(X)
  n <- nrow(Y)
  J <- ncol(Y)
  
  # get subset 
  if (is.null(reference_set)) {
    ref_set <- try(get("reference_set", envir = environment(constraint_fn[[k]])))
    if (inherits(ref_set, "try-error")) {
      stop()
    } else {
      ref_set_k <- ref_set[[k_constr]]
    }
  } else {
    ref_set_k <- reference_set
  }
  other_set <- (1:J)[-ref_set_k]
  
  if (is.list(constraint_fn)) {
    constraint_fn <- constraint_fn[[k_constr]]
    constraint_grad_fn <- constraint_grad_fn[[k_constr]]
  }
  
  # data frame to hold results
  it_df <- data.frame(it = 1:maxit,
                      lik = NA,
                      test_stat = NA)
  
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
  
  # do this only for taxa in subset 

  above_g <- setdiff(
    which(B[k_constr, ] > constraint_fn(B[k_constr, ]) & 
            (1:J) %in% ref_set_k),
    c(j_ref, j_constr)
  )
  above_g <- above_g[order(B[k_constr, above_g])]
  below_g <- setdiff(
    which(B[k_constr, ] < constraint_fn(B[k_constr, ]) &
            (1:J) %in% ref_set_k),
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
  if (verbose) {
    print(paste0("Starting ll = ", ll))
  }
  test_stat <- NA

  #initialize 'keep_going' indicator for continuing optimization as well as 'iter'
  #to track optimization iterations
  keep_going <- TRUE
  iter <- 1

  # temporarily set use_optim to `TRUE` always
  # use_optim <- TRUE

  lik_change <- rep(NA, null_window - 1)
  test_stat_prop_change <- rep(NA, null_window - 1) 
  
  #iterate until something happens
  while (keep_going) {
    prev_B <- B
    prev_ll <- ll
    prev_test_stat <- test_stat

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
            scale <- max(abs(Matrix::diag(info)))
            repeat {
              if (lambda > 0) {
                info_reg <- info +
                  lambda *
                    #Matrix::Diagonal(
                    #  n = nrow(info),
                    #  x = pmax(abs(Matrix::diag(info)), 1)
                    #)
                  Matrix::Diagonal(n = nrow(info), x = scale)
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
            
            # add in max step 
            max_step_dir <- max(abs(step_dir), na.rm = TRUE)
            if(max_step_dir > max_step){
              step_dir <- step_dir / max_step_dir * max_step
            }
            
            # -- back-tracking Armijo line search
            step <- 1
            f_curr <- fn(theta)
            slope <- sum(grad * step_dir)
            
            for (bt in seq_len(backtrack_max)) {
              theta_new <- theta + step * step_dir
              f_new <- fn(theta_new)
              
              # added in check for if f_new is getting infinite
              if (is.finite(f_new) && f_new <= f_curr + c1 * step * slope) {
                break
              }
              step <- step / 2
            }
            
            # if f_new is always infinite, set step to 0 to hit convergence, move on without updating
            if (!is.finite(f_new)) {
              step <- 0
              f_new <- f_curr
              theta_new <- theta
            }

            # -- update
            theta <- theta_new
            
            # if (verbose) {
            #   cat(sprintf(
            #     "Iter %2d  f = %-12.6g  abs theta change = %-10.3g  lambda = %g\n",
            #     iter,
            #     f_new,
            #     max(abs(step * step_dir)),
            #     lambda
            #   ))
            # }

            # -- convergence?
            if (max(abs(step * step_dir)) < tol) {
              
              if (step == 0) {
                return(list(
                  par = theta,
                  value = f_new,
                  converged = FALSE,
                  iter = iter
                ))
              }
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
    
    # loop through taxa not in the reference set 
    for (j_left in other_set) {
      
      # for taxa not in reference set, update with same approach as in emuFit_micro
      update <- micro_fisher(
        X,
        Yj = Y[, j_left, drop = FALSE],
        Bj = B[, j_left, drop = FALSE],
        z,
        stepsize = max_step,
        c1 = c1
      )
      
      B[, j_left] <- B[, j_left] + update
      
      #update z
      z <- update_z(Y = Y, X = X, B = B)
      
    }

    #compute ll again
    log_means <- X %*% B
    for (i in 1:nrow(Y)) {
      log_means[i, ] <- log_means[i, ] + z[i]
    }

    ll <- sum(Y * log_means - exp(log_means))
    score_res <- try(get_score_stat(X = X, Y = Y, X_cup = X_cup, B = B, k_constr = k_constr,
                                    j_constr = j_constr, constraint_grad_fn = constraint_grad_fn,
                                    indexes_to_remove = (j_ref - 1)*p + 1:p, j_ref = j_ref, J = J,
                                    n = n, p = p))
    if (inherits(score_res, "try-error")) {
      stop(paste0("fit_null_symmetric() failed because the test statistic could not be computed for iteration ", iter))
    } else {
      test_stat <- score_res$score_stat
    }

    it_df$lik[iter] <- ll
    it_df$test_stat[iter] <- test_stat

    if (iter > 1) {
      test_stat_change <- abs((test_stat - prev_test_stat)/prev_test_stat)
    } 
    ll_change <- abs(ll - prev_ll) / abs(prev_ll)
    
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
      message("ll increased by ", round(ll_change * 100, 2), "%")
      if (iter > 1) {
        message("test statistic changed by ", round(test_stat_change * 100, 2), "%")
      }
      #message("gr_norm = ", round(sum(gr^2), 1))
      #message("max abs diff in B = ", signif(max(abs(B - prev_B)), 2))
    }
    
    if (iter %in% 2:null_window) {
      lik_change[iter-1] <- ll_change
      test_stat_prop_change[iter-1] <- test_stat_change
    } else if (iter > null_window) {
      lik_change[1:(null_window - 2)] <- lik_change[2:(null_window - 1)]
      lik_change[null_window - 1] <- ll_change
      test_stat_prop_change[1:(null_window - 2)] <- test_stat_prop_change[2:(null_window - 1)]
      test_stat_prop_change[null_window - 1] <- test_stat_change 
    }

    #recompute above_g and below_g
    above_g <- setdiff(
      which(B[k_constr, ] > constraint_fn(B[k_constr, ]) &
              (1:J) %in% ref_set_k),
      c(j_ref, j_constr)
    )
    above_g <- above_g[order(B[k_constr, above_g])]
    below_g <- setdiff(
      which(B[k_constr, ] < constraint_fn(B[k_constr, ]) &
              (1:J) %in% ref_set_k),
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

    # update stopping criteria
    # if (max(abs(B - prev_B)) < B_tol) {
    #   keep_going <- FALSE
    #   converged <- TRUE
    # }
    if (!(ignore_stop) & (iter - 1) >= null_window) {
      if (max(lik_change) < tol_lik & max(test_stat_prop_change) < tol_test_stat) {
        keep_going <- FALSE
        converged <- TRUE
      }
    }
  }
  
  it_df <- it_df[1:(iter - 1), ]

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
  
  # check if converged if ignoring stopping criteria
  if (ignore_stop) {
    if (max(lik_change) < tol_lik & max(test_stat_prop_change) < tol_test_stat) {
      converged <- TRUE
    }
  }
  
  #return a list containing B 'n' fixin's
  return(list(
    B = B,
    k_constr = k_constr,
    j_constr = j_constr,
    niter = iter,
    converged = converged,
    Bs = df,
    it_df = it_df
  ))
}
