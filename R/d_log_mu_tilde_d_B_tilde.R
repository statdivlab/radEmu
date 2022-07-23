

d_log_mu_tilde_d_B_tilde <- function(Y_tilde, X_tilde, log_mu_tilde,n,J){


  Us <- lapply(1:n,
               function(i) diag(rep(1,J)) -
                 matrix(1,nrow = J, ncol = 1)%*%matrix(exp(log_mu_tilde[(i - 1)*J + 1:J,]),nrow = 1)
  )
  if(n >1){
  U <- bdiag_m(Us)
  } else{
    U <- Us[[1]]
  }

  return(U%*%X_tilde)
}
