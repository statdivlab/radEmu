

eta_jacobian <- function(beta_tilde_J,
                         tau,
                         X_tilde_J,
                         J,
                         i,
                         n){
  jacobian <- matrix(0,nrow = J,ncol = ncol(X_tilde_J) + n)

    which_indices <- 1:J + (i - 1)*J
    p_i <- exp(X_tilde_J[which_indices,]%*%beta_tilde_J)
    p_i <- p_i/sum(p_i)
    premult_mat <- matrix(1,nrow = J, ncol = 1)%*%Matrix::t(p_i)
    premult_mat <- diag(rep(1,J)) - premult_mat
    d_eta_i_d_beta_tilde <- premult_mat%*%X_tilde_J[which_indices,]
    e_i <- matrix(0,nrow = 1, ncol = n)
    e_i[1,i] <- 1
    d_eta_i_d_tau <- (1/tau[i])*matrix(1,nrow= J, ncol= 1)%*%e_i

    jacobian <- jacobian + cbind(d_eta_i_d_beta_tilde, d_eta_i_d_tau)

  return(jacobian)
}
