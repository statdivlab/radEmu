

asymmetric_hat <- function(beta_tilde_J,
                           tau,
                           X_tilde_J,
                           J,
                           n,
                           p,
                           collect_jacobians = TRUE,
                           verbose= FALSE){

  expected_info_list <- expected_info_repar(beta_tilde_J,
                                       tau,
                                       X_tilde_J,
                                       J,
                                       n,
                                       p,
                                       verbose = verbose)
  expected_info <- expected_info_list$expected_info
  # expected_info_eigen <- eigen(expected_info)

  expected_info_eigen <- eigen(expected_info[1:(p*(J - 1)),1:(p*(J - 1))])


  expected_info_half_inv <-  expected_info_eigen$vectors %*%
    diag(sapply( expected_info_eigen$values,function(x) ifelse(x>0,1/sqrt(x),0))) %*%
    t(expected_info_eigen$vectors)

  tau_start <- p*(J - 1) + 1
  tau_end <- nrow(expected_info)
  expected_info_half_inv <- Matrix::bdiag(expected_info_half_inv,
                                          Matrix::Diagonal(
                                            x = 1/sqrt(diag(expected_info)[tau_start:tau_end])))

  augmentations <- matrix(0,nrow = n, ncol = J)

  for(i in 1:n){
    if(verbose){
      message("Computing data augmentation for observation ",i,".")
    }
    which_indices <- 1:J + (i - 1)*J
    jacobian_i <- expected_info_list$jacobians[[i]]
    mu_i <- expected_info_list$mus[[i]]
    postmult_jacobian_i <- jacobian_i%*%  expected_info_half_inv
    H_i <- Matrix::tcrossprod(postmult_jacobian_i)%*%diag(as.numeric(mu_i))
    augmentations[i,] <- diag(as.matrix(H_i))[1:J]/2
  }

  return(augmentations)


}
