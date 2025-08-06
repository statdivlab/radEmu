compare_null <- function(maxit = 1000, record_gap = 10, X, Y, 
                         k_constr, j_constr, verbose = F, alt_fit = NULL) {
  
  # hyperparams
  n <- nrow(X)
  J <- ncol(Y)
  p <- ncol(X)
  
  # set constraint function and its gradient 
  constraint_fn <- function(x) {pseudohuber_median(x, 0.1)}
  constraint_grad_fn <- function(x) {dpseudohuber_median_dx(x, 0.1)}
  
  X_cup <- X_cup_from_X(X, J)
  j_ref <- get_j_ref(Y)
  
  # if alt fit is NULL, get this
  if (is.null(alt_fit)) {
    alt_fit <- emuFit(X = X, Y = Y, compute_cis = FALSE, run_score_tests = FALSE)
    if (verbose) {
      message("Fitting model under alternative")
    }
  }
  Y <- alt_fit$Y_augmented
  
  # make result data frames for old and new algorithms 
  record_its <- (1:maxit)[1:maxit %% record_gap == 0]
  old <- data.frame(it = record_its,
                    time = NA,
                    log_lik = NA,
                    mse_lag = NA,
                    max_abs_B = NA,
                    constraint_diff = NA,
                    test_stat = NA)
  new <- old 
  new$alg <- "new"
  old$alg <- "old"
  
  # fit old null
  message("Fitting old verison of null fit")
  start_old <- proc.time()
  
  # copying and modifying from fit_null (defaults from emuFit): 
  rho_init = 1; tau = 2; kappa = 0.8; inner_tol = 1; constraint_tol = 1e-5
  max_step = 1; c1 = 1e-4; inner_maxit = 25
  
  B <- alt_fit$B
  
  #store cols of B to update in a vector
  j_to_update <- (1:J)[(1:J) != j_ref]
  #in vector format, which indexes to remove (bc B^j_constr = 0 id. constr.)
  indexes_to_remove <- (j_ref - 1)*p + 1:p
  
  
  #change id. constr. by subtracting approp. col of B from others
  for(k in 1:p) {
    B[k,] <- B[k,] - B[k,j_ref]
  }
  
  #update z
  z <- update_z(X = X, Y = Y, B = B)
  
  log_mean <- X%*%B +
    matrix(z,ncol = 1)%*%matrix(1,ncol = J, nrow = 1)
  
  Bs <- NULL
    
  #set iteration to zero
  iter <- 0
  
  #compute gap (i.e. g(B_k) - B_kj)
  gap <- constraint_fn(B[k_constr,]) - B[k_constr,j_constr]
  init_gap <- abs(gap)
  
  #set rho equal to initial value
  rho <- rho_init
  
  #initiate u
  u <- rho*gap
  
  loop_j_to_update <- j_to_update
  
  use_max_inner_it <- FALSE
  B_diff <- Inf
  
  while(iter <= maxit & !is.infinite(rho)) {
    
    iter <- iter + 1    #increment iteration
    inner_iter <- 0     #initiate internal iteration
    
    #evaluate augmented Lagrangian
    log_means <- do.call(cbind,lapply(1:J,function(j) X%*%B[,j] + z))
    
    #  get current value of augmented lagrangian
    curr_lag_val <- sum(-log_means*Y + exp(log_means)) + u*gap + (rho/2)*(gap^2)
    
    old_gap <- gap # old value of g(B_k) - B_kj
    old_B <- B #old value of B
    inner_diff <- Inf #initiate "observed" max abs val diff in B between inner iterations
    #at Inf so inner loop runs at least once
    #inner loop:
    while(((inner_diff > inner_tol | use_max_inner_it) & inner_iter <= inner_maxit)| inner_iter == 0) {
      
      inner_old_B <- B
      
      #perform block update to B as described in null estimation section of Clausen & Willis (2024)
      update <- macro_fisher_null(X = X,
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
                                  c1 = 1e-4)
      B <- B + update$update
      gap <- update$gap
      z <- update_z(Y,X,B)
      
      #did we move much?
      inner_diff <- max(abs(B - inner_old_B))
      
      #increment inner iteration counter
      inner_iter <- inner_iter + 1
      
    }
    #did we move much since last *outer* loop iteration?
    B_diff <- max(abs(B - old_B))
    
    #update u and rho
    if ( abs(gap) > constraint_tol) {
      u <- u + rho*gap
      if (abs(gap/old_gap) > kappa) {
        rho <- rho*tau
      } 
    }
    
    if (abs(gap) < constraint_tol) {
      use_max_inner_it <- TRUE
    } else{
      use_max_inner_it <- FALSE
    }
    
    if (iter %% record_gap == 0) {
      if (verbose) {
        message(paste0("old alg, iteration ", iter))
      }
      ind <- which(old$it == iter)
      old$time[ind] <- (proc.time() - start_old)[3]
      old$max_abs_B[ind] <- B_diff
      old$constraint_diff[ind] <- abs(gap)
      
      log_means <- X%*%B
      for(i in 1:n){
        log_means[i,] <- log_means[i,] + z[i]
      }
      old$log_lik[ind] <- sum(Y * log_means - exp(log_means))
      old$mse_lag[ind] <- min_mse_lag(X = X, Y = Y,
                                      B = B, constraint_grad_fn = constraint_grad_fn,
                                      k_constr = k_constr, j_constr = j_constr, j_ref = j_ref)
      
      old$test_stat[ind] <- get_score_stat(X = X, Y = Y, X_cup = X_cup, B = B, k_constr = k_constr,
                                           j_constr = j_constr, constraint_grad_fn = constraint_grad_fn,
                                           indexes_to_remove = indexes_to_remove, j_ref = j_ref, J = J,
                                           n = n, p = p)$score_stat
    }
    
  }
  
  # fit new null
  message("Fitting new version of null fit")
  start_new <- proc.time() 
  
  # copied and modified from fit_null_symmetric:
  B <- alt_fit$B
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
      
      # Fisher scoring (by Amy & copilot)
      
      #------------------------------------------------------
      # 1.  Build problem-specific closures once
      # -----------------------------------------------------
      make_fs_problem <- function(js,
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
      fisher_scoring <- function(theta0,
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
      prob <- make_fs_problem(js,
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
      
      
      #update z
      z <- update_z(Y = Y, X = X, B = B)
    }
    
    #compute ll again
    log_means <- X %*% B
    for (i in 1:nrow(Y)) {
      log_means[i, ] <- log_means[i, ] + z[i]
    }
    
    ll <- sum(Y * log_means - exp(log_means))
    
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
    
    # save relevant pieces
    if (iter %% record_gap == 0) {
      if (verbose) {
        message(paste0("new alg, iteration ", iter))
      }
      ind <- which(new$it == iter)
      new$time[ind] <- (proc.time() - start_new)[3]
      new$max_abs_B[ind] <- max(abs(B - prev_B))
      new$log_lik[ind] <- ll
    
      new$mse_lag[ind] <- min_mse_lag(X = X, Y = Y,
                                      B = B, constraint_grad_fn = constraint_grad_fn,
                                      k_constr = k_constr, j_constr = j_constr, j_ref = j_ref)
      
      new$test_stat[ind] <- get_score_stat(X = X, Y = Y, X_cup = X_cup, B = B, k_constr = k_constr,
                                           j_constr = j_constr, constraint_grad_fn = constraint_grad_fn,
                                           indexes_to_remove = indexes_to_remove, j_ref = j_ref, J = J,
                                           n = n, p = p)$score_stat
    }
    
    #increment iteration counter
    iter <- iter + 1
    
    #determine whether we should stop iterating
    if (iter > maxit) {
      keep_going <- FALSE
      converged <- FALSE
    }
    
  }
  
  #if we changed j_ref, change it back
  if (!is.null(arch_j_ref)) {
    for (k in 1:p) {
      B[k, ] <- B[k, ] - B[k, arch_j_ref]
    }
  }
  
  # return results
  res <- rbind(new, old)
  return(res)
  
}