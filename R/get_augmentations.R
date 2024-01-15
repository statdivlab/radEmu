

get_augmentations <- function(X,
                              G,
                              Y,
                              B){


  #dimensions
  n <- nrow(Y)
  p <- nrow(B)
  J <- ncol(Y)

  #parametrize so last column of B is dropped
  for(k in 1:p){
    B[k, ] <- B[k, ] - B[k, J]
  }
  B_tilde_J <- B_cup_from_B(B)
  to_delete <- p*(J - 1) + 1:p
  B_tilde_J <- B_tilde_J[-to_delete, , drop = FALSE]


  Y_rowsums <- rowSums(Y)
  z_tilde <- numeric(n)
  for(i in 1:n){
    z_tilde[i] <- log(Y_rowsums[i]) - sum_of_logs(X[i, , drop= FALSE] %*% B)
  }

  theta_tilde <- rbind(B_tilde_J, matrix(z_tilde, ncol = 1))

  log_means <- G %*% theta_tilde

  W <- Matrix::Diagonal(x = exp(log_means@x))

  info <- Matrix::crossprod(G, W) %*% G
  info <- methods::as(info, "symmetricMatrix")
  # info <- as(info,"dsTMatrix")
  # info <- spam::as.spam.dgCMatrix(info)
  info_chol <- Matrix::chol(info, pivot = FALSE)
  info_half_inv <- Matrix::solve(info_chol)
  # info_half_inv <- info_chol_inv
  # info_half_inv <- 0.5*(info_chol_inv + t(info_chol_inv))
  # eigen(info_chol)
  # info_schur <- Matrix::Schur(info)
  # info_eigen <- eigen(info, symmetric = TRUE)
  # info_eigen <- spam::eigen.spam(info,
  #                                symmetric = TRUE,
  #                                nev = nrow(info))
  # info_half_inv <-  info_eigen$vectors  %*%
  #   Matrix::tcrossprod(Matrix::Diagonal(x = 1/sqrt(info_eigen$values)), info_eigen$vectors)
  # info_half_inv <- Matrix::t(info_chol_inv)
  augmentations <- matrix(0, nrow = n, ncol =J)
  # pv <- profvis::profvis({
  for(i in 1:n){
    # print(i)
    G_i <- G[1:J + (i - 1)*J, , drop = FALSE]
    G_i <- G_i %*% info_half_inv

    # G_i %*%
    # for(j in 1:J){
    #   augmentations[i, j] <- sum(G_i[j, ]^2)*exp(log_means[(i - 1)*J +j, ])/2
      augmentations[i, ] <-
      Matrix::rowSums(Matrix::Diagonal(x = exp(log_means[(i - 1)*J + 1:J, ])) %*% (G_i^2) )/2
  #     G_i
  #   }
  }
  # })


  return(augmentations)



}
