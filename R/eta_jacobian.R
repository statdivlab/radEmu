

eta_jacobian <- function(beta_tilde_J,
                         tau,
                         X_tilde_J,
                         J,
                         i,
                         n){

  which_indices <- 1:J + (i - 1)*J
  p_i <- exp(X_tilde_J[which_indices,]%*%beta_tilde_J)
  p_i <- p_i/sum(p_i)
  premult_mat <- Matrix::Matrix(1,nrow = J, ncol = 1)%*%Matrix::t(p_i)
  premult_mat <- Matrix::Diagonal(x = rep(1,J)) - premult_mat
  d_eta_i_d_beta_tilde <- premult_mat%*%X_tilde_J[which_indices,]
  # d_eta_i_d_beta_tilde <- as(d_eta_i_d_beta_tilde,"sparseMatrix")
  # e_i <- Matrix(0,nrow = 1, ncol = n)
  # e_i[1,i] <- 1
  # d_eta_i_d_tau <- Matrix::Matrix(0,nrow = J, ncol = n)
  # d_eta_i_d_tau[,i] <- 1/tau[i]

  d_eta_i_d_tau <- Matrix::sparseMatrix(i = 1:J,
                                        j = rep(i, J),
                                        x = rep(1/tau[i],J),
                                        dims = c(J, n))

  jacobian <- cbind(d_eta_i_d_beta_tilde, d_eta_i_d_tau)
  # jacobian <- list()
  # jacobian <- list("dbeta" = d_eta_i_d_beta_tilde,
  #                  "dtau" = d_eta_i_d_tau)



  return(jacobian)
}
