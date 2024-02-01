
#perform score test
score_test <- function(B, #B (MPLE)
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
                       return_both_score_pvals = FALSE){
  n <- nrow(Y)
  J <- ncol(Y)
  p <- ncol(X)


  tries_so_far <- 0
  accept_try <- FALSE
  good_enough_fit <- FALSE
  while(!accept_try){
    #fit under null
  constrained_fit <- try(fit_null(B = B, #B (MPLE)
                                      Y = Y, #Y (with augmentations)
                                      X = X, #design matrix
                                      X_cup = X_cup,
                                      k_constr = k_constr, #row index of B to constrain
                                      j_constr = j_constr, #col index of B to constrain
                                      constraint_fn = constraint_fn, #constraint function
                                      constraint_grad_fn = constraint_grad_fn, #gradient of constraint fn
                                      # constraint_hess_fn = constraint_hess_fn,
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
                                      # I = I,
                                      # Dy = Dy
  ))
  if(inherits(constrained_fit,"try-error")){
    accept_try <- FALSE
  } else{
    if((abs(constrained_fit$gap) <= constraint_tol) &
       (constrained_fit$niter < maxit)){
      accept_try <- TRUE
      good_enough_fit <- TRUE
    } else{
      tau <- tau^(3/4)
      inner_maxit <- 2*inner_maxit
      message("Constrained optimization failed to converge within iteration limit;
retrying with smaller penalty scaling parameter tau and larger inner_maxit.")
    }
  }
  tries_so_far <- tries_so_far + 1
  if(tries_so_far == ntries){
    accept_try <- TRUE
  }
  }


  if(!good_enough_fit){
    warning("Optimization for null fit with k = ",k_constr," and j = ",j_constr," failed to converge across ", ntries, ifelse(ntries>1," attempts."," attempt."))
}
  B <- constrained_fit$B
  z <- update_z(Y,X,B)
  p <- ncol(X)
  J <- ncol(Y)
  n <- nrow(Y)

  # message("Make X_cup an argument of score_test and calculate it once in emuFit")
  # X_cup = X_cup_from_X(X,J)
    #indexes in long format corresponding to the j_constr-th col of B
  #get score stat
  indexes_to_remove <- (j_ref - 1)*p + 1:p
  score_stat <-
    get_score_stat(Y = Y,
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
                 Dy = Dy)


  score_stat <- score_stat*(n/(n - 1))

 
  if(!return_both_score_pvals){ #typically we want only one score p-value
                                #(using only one version of information matrix)

  return(list("score_stat" = score_stat,
              "pval" = pchisq(score_stat,1,lower.tail = FALSE),
              "log_pval" = pchisq(score_stat,1,lower.tail = FALSE, log.p = TRUE),
              "niter" = constrained_fit$niter,
              "convergence" = ifelse(constrained_fit$niter>=maxit,'iteration limit reached','converged'),
              "proj_score" = constrained_fit$proj_score,
              "gap" = constrained_fit$gap,
              "u" = constrained_fit$u,
              "rho" = constrained_fit$rho,
              "null_B" = constrained_fit$B,
              # "score_stats" = constrained_fit$score_stats,
              "Bs" = constrained_fit$Bs))
  } else{ 
    #for simulations -- if we want to return both the score p-value using
    #information from full model fit and from null model
      score_stat_with_null_info <-
        get_score_stat(Y = Y,
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
                       Dy = Dy)

      score_stat_with_null_info <- score_stat_with_null_info*(n/(n - 1))

      return(list("score_stat" = score_stat,
                  "pval" = pchisq(score_stat,1,lower.tail = FALSE),
                  "log_pval" = pchisq(score_stat,1,lower.tail = FALSE, log.p = TRUE),
                  "score_stat_null_info" = score_stat_with_null_info,
                  "pval_null_info" = pchisq(score_stat_with_null_info,1,lower.tail = FALSE),
                  "log_pval_null_info" = pchisq(score_stat_with_null_info,1,lower.tail = FALSE,log.p = TRUE),
                  "niter" = constrained_fit$niter,
                  "convergence" = ifelse(constrained_fit$niter>=maxit,'iteration limit reached','converged'),
                  "proj_score" = constrained_fit$proj_score,
                  "gap" = constrained_fit$gap,
                  "u" = constrained_fit$u,
                  "rho" = constrained_fit$rho,
                  "null_B" = constrained_fit$B,
                  # "score_stats" = constrained_fit$score_stats,
                  "Bs" = constrained_fit$Bs))
    }


}
