
#get design matrix in beta_vec (w.o. B^J) and z
get_G_for_augmentations <- function(X,
                                    J,
                                    n,
                                    X_cup){
  p <- ncol(X)
  to_delete <- p*(J - 1) + 1:p
  # X_tilde_J <- X_cup_from_X(X,J)
  X_tilde_J <- X_cup
  for(i in 1:n){
    X_tilde_J[(i - 1)*J  + J,] <- 0
  }

  X_tilde_J <- X_tilde_J[,-to_delete]

  Z <- Matrix::sparseMatrix(i = 1:(n*J),
                            j = rep(1:n, each = J),
                            x = rep(1,n*J))
  G <- cbind(X_tilde_J,Z)
}
