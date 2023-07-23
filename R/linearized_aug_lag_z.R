
linearized_aug_lag_z <-
  function(X,
           Y,
           j,
           j_constr,
           k_constr,
           B,
           z,
           B0,
           u,
           rho,
           constraint_fn_at_B0,
           constraint_grad_at_B0,
           constraint_fn,
           constraint_grad_fn,
           compute_gradient = FALSE,
           compute_hessian = FALSE){

    lin_gap <- B[k_constr,j_constr] -
      constraint_fn(B[k_constr,])

    log_means <- X%*%B[,j,drop = FALSE] + z
    lin_lag_val <-
      sum(-Y[,j,drop = FALSE]*log_means +
                       exp(log_means)) +
      u*lin_gap +
      (rho/2)*(lin_gap)^2
    if(compute_gradient){
      gr_val <- as.numeric(-t(X)%*%(Y[,j,drop = FALSE] - exp(log_means)))
      gr_val[k_constr] <-  gr_val[k_constr] +
        (u + rho*lin_gap)*(as.numeric(j == j_constr) -    constraint_grad_fn(B[k_constr,])[j])
    } else{
      gr_val <- NULL
    }

    if(compute_hessian){
      hess_val <- t(X)%*%diag(as.numeric(exp(log_means)))%*%X
      hess_val[k_constr,k_constr] <-     hess_val[k_constr,k_constr] +
        rho*(ifelse(j == j_constr,1,0) - constraint_grad_fn(B[k_constr,])[j])^2

    } else{
      hess_val <- NULL
    }

    return(list("value" = lin_lag_val,
                "gr"= gr_val,
                "hess" = hess_val))
  }
