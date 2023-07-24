

micro_score_test <- function(Y,
                             X,
                             B = NULL,
                             constraint_fn,
                             constraint_grad_fn,
                             null_k,
                             null_j,
                             maxit = 500,
                             tolerance = 1e-2,
                             verbose = TRUE,
                             step_ratio = 0.5,
                             rho_init = 1,
                             rho_scaling = 2,
                             gap_tolerance = 1e-3){
  n <- nrow(Y)
  J <- ncol(Y)
  p <- ncol(X)

  # constrained_fit <- emuFit_micro_constrained(X = X,
  #                          Y = Y,
  #                          k_constr = null_k,
  #                          j_constr = null_j,
  #                          B=  B,
  #                          maxit = maxit,
  #                          constraint_fn = constraint_fn,
  #                          constraint_grad_fn = constraint_grad_fn,
  #                          tolerance = tolerance,
  #                          verbose = verbose
  #                          )

  constrained_fit <-
    lagrangian_fit(X = X,
                   Y = Y,
                   B = B,
                   constraint_fn = constraint_fn,
                   constraint_grad_fn = constraint_grad_fn,
                   j_constr = null_j,
                   k_constr = null_k,
                   maxit = 250,
                   rho_init = rho_init,
                   rho_scaling = rho_scaling,
                   tolerance = tolerance,
                   gap_tolerance = gap_tolerance,
                   c1 = 1e-4,
                   verbose = verbose)

  B <- constrained_fit
  z <- update_z_no_wts(Y,X,B)

  X_cup = X_cup_from_X(X,J)
  scores <- vector(n,mode = "list")

  #compute score contributions of observations i = 1 through n
 for(i in 1:n){
   X_cup_i <- X_cup[(i - 1)*J + 1:J,]
   scores[[i]] <- as.matrix(dpll_dB_cup(X[i,,drop = FALSE],Y[i,,drop = FALSE],B))
 }

    H <- matrix(0,nrow = p, ncol = J)
    for(j in 1:J){
      # H[null_k,j] <- get_dg_dBj(B = B,
      #                           j = j,
      #                           k_star = null_k,
      #                           huber_param = huber_param,
      #                           constraint_type = constraint_type,
      #                           for_testing = TRUE)[null_k]

      H[null_k,j] <- constraint_grad_fn(B[null_k,])[j]
    }
    H[null_k,null_j] <- H[null_k,null_j] -1
    H_cup <- B_cup_from_B(H)

    B_cup <- B_cup_from_B(B)
    I <- f_info(Y,B_cup,B,X,X_cup)
    I_geninv <- eigen(I)
    n_eigen <- p*(J - 1)
    n_remaining <- p*J - n_eigen
    I_geninv <- I_geninv$vectors%*%diag(c(1/I_geninv$values[1:n_eigen],rep(0,n_remaining)))%*%t(I_geninv$vectors)
    H_cup <- as.matrix(H_cup)
    H_I_inv <- t(H_cup)%*%I_geninv

    Umat <- t(do.call(cbind,scores))
    Dy <- Reduce("+",lapply(scores,function(x) tcrossprod(x)))
    score <- Reduce("+",scores)
    score <- as.matrix(score)

    outside <- H_I_inv%*%score

    inside <- H_I_inv%*%Dy%*%t(H_I_inv)

    score_stat <- (n/(n - 1))*
      as.numeric(outside)^2/as.numeric(inside)

    return(list("score_stat" = score_stat, pval = pchisq(score_stat,1,lower.tail = FALSE)))


}
