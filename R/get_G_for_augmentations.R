
#get design matrix in beta_vec (w.o. B^J) and z
get_G_for_augmentations <- function(X,
                                    J,
                                    n,
                                    X_cup){
  p <- ncol(X)
  #calculate indices of columns to remove from X_cup (aka X_tilde, the expanded design matrix)
  to_delete <- p*(J - 1) + 1:p #col indices for elements of B^J in long B format
  X_tilde_J <- X_cup
  for(i in 1:n){
    X_tilde_J[(i - 1)*J  + J,] <- 0
  }

  X_tilde_J <- X_tilde_J[,-to_delete]

  #design matrix in z
  Z <- Matrix::sparseMatrix(i = 1:(n*J),
                            j = rep(1:n, each = J),
                            x = rep(1,n*J))
  G <- cbind(X_tilde_J,Z)
  return(G)
}
