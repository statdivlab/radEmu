
#get data augmentations for Firth penalized estimation
#general method described in Kosmidis and Firth (2011)
#specifics here described in estimation section of Clausen and Willis (2024)
get_augmentations <- function(X,
                              G, #design matrix in beta_vec and z
                              Y,
                              B){


  #dimensions
  n <- nrow(Y)
  p <- nrow(B)
  J <- ncol(Y)

  #parametrize so last column of B is / can be dropped
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

  #information matrix (we're using tricky equivalency to poisson model)
  info <- Matrix::crossprod(G, W) %*% G
  # info <- methods::as(info, "symmetricMatrix")
  info <- Matrix::forceSymmetric(info)
  #cholesky decomposition
  info_chol <- Matrix::chol(info, pivot = FALSE)
  #info^(-0.5) by inverting cholesky decomp
  info_half_inv <- Matrix::solve(info_chol)
 
  augmentations <- matrix(0, nrow = n, ncol =J)

  #get augmentations
  for(i in 1:n){
    # print(i)
    G_i <- G[1:J + (i - 1)*J, , drop = FALSE]
    G_i <- G_i %*% info_half_inv

      augmentations[i, ] <-
      Matrix::rowSums(Matrix::Diagonal(x = exp(log_means[(i - 1)*J + 1:J, ])) %*% (G_i^2) )/2

  }



  return(augmentations)



}

# this actually isn't faster in tests
get_augmentations_par <- function(X,
                                  G, #design matrix in beta_vec and z
                                  Y,
                                  B,
                                  par = FALSE){
  
  
  #dimensions
  n <- nrow(Y)
  p <- nrow(B)
  J <- ncol(Y)
  
  #parametrize so last column of B is / can be dropped
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
  
  #information matrix (we're using tricky equivalency to poisson model)
  info <- Matrix::crossprod(G, W) %*% G
  # info <- methods::as(info, "symmetricMatrix")
  info <- Matrix::forceSymmetric(info)
  #cholesky decomposition
  info_chol <- Matrix::chol(info, pivot = FALSE)
  #info^(-0.5) by inverting cholesky decomp
  info_half_inv <- Matrix::solve(info_chol)
  
  
  
  if (.Platform$OS.type == "unix" & par) {
    
    cores <- parallel::detectCores()
    aug_res <- parallel::mclapply(1:n, function(i) {
      G_i <- G[1:J + (i - 1)*J, , drop = FALSE]
      G_i <- G_i %*% info_half_inv
      Matrix::rowSums(Matrix::Diagonal(x = exp(log_means[(i - 1)*J + 1:J, ])) %*% (G_i^2) )/2
    }, mc.cores = cores - 1)
    augmentations <- do.call(rbind, aug_res)
    
  } else {
    augmentations <- matrix(0, nrow = n, ncol = J)
    for(i in 1:n){
      G_i <- G[1:J + (i - 1)*J, , drop = FALSE]
      G_i <- G_i %*% info_half_inv
      
      augmentations[i, ] <-
        Matrix::rowSums(Matrix::Diagonal(x = exp(log_means[(i - 1)*J + 1:J, ])) %*% (G_i^2) )/2
      
    }
  }
  
  return(augmentations)
  
}