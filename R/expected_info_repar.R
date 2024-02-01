
#another function to get expected information (and possibly jacobians of 
#mean function in parameters)
expected_info_repar <- function(beta_tilde_J,
                                tau,
                                X_tilde_J,
                                J,
                                n,
                                p,
                                collect_byproducts = TRUE,
                                verbose = FALSE){

  expected_info <- Matrix::Matrix(0,nrow = n + p*(J - 1), ncol = n + p*(J - 1))

  if(collect_byproducts){
    jacobians <- vector(n, mode = "list")
    mus <- vector(n,mode = "list")
  } else{
    jacobians <- NULL
    mus <- NULL
  }
  # info_bits <- vector(n,mode = "list")

  # pv5 <- profvis({
  for(i in 1:n){
    if(verbose){
      message("Computing information contribution of observation ", i,".")
    }

    which_indices <- 1:J + (i - 1)*J
    d_eta_i_d_theta_J <- eta_jacobian(beta_tilde_J,
                                      tau,
                                      X_tilde_J,
                                      J,
                                      i,
                                      n)
    if(collect_byproducts){
      jacobians[[i]] <- d_eta_i_d_theta_J
    }
    mu_i <- exp(X_tilde_J[which_indices,]%*%beta_tilde_J + log(tau[i]) -
                  log(sum(exp(X_tilde_J[which_indices,]%*%beta_tilde_J ))))
    mus[[i]] <- mu_i
    expected_info <- expected_info +
      Matrix::crossprod(d_eta_i_d_theta_J,Matrix::crossprod(Matrix::Diagonal(x = as.numeric(mu_i)),d_eta_i_d_theta_J))
    # Matrix::t(d_eta_i_d_theta_J)%*%diag(as.numeric(mu_i)) %*%
    # d_eta_i_d_theta_J
  }
  # })
  return(list("expected_info" = expected_info,
              "jacobians" = jacobians,
              "mus" = mus))
}
