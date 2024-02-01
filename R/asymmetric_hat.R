
#this function creates the "asymmetric hat" matrix used for 
#creating data augmentations as described in Kosmidis & Firth (2011)

asymmetric_hat <- function(beta_tilde_J, #beta expressed as a vector
                           tau, #column sums of Y (without data augmentations)
                           X_tilde_J, #expanded design matrix (for Y as a single vector)
                           J, #number of taxa
                           n, #number of samples
                           p, #number of predictors, incl. intercept
                           collect_jacobians = TRUE,
                           verbose= FALSE){

  #get the expected information matrix 
  expected_info_list <- expected_info_repar(beta_tilde_J,
                                       tau,
                                       X_tilde_J,
                                       J,
                                       n,
                                       p,
                                       verbose = verbose)
  expected_info <- expected_info_list$expected_info
  # expected_info_eigen <- eigen(expected_info)

  #eigendecompose expected information after removing J-th taxon's parameters 
  #(implicitly imposing B^J = 0)
  expected_info_eigen <- eigen(expected_info[1:(p*(J - 1)),1:(p*(J - 1))])

  # get info^(-0.5) using eigendecomposition
  expected_info_half_inv <-  expected_info_eigen$vectors %*%
    diag(sapply( expected_info_eigen$values,function(x) ifelse(x>0,1/sqrt(x),0))) %*%
    t(expected_info_eigen$vectors)

  # now we add on information in tau
  tau_start <- p*(J - 1) + 1 #index where tau starts 
  tau_end <- nrow(expected_info) #index where tau ends
  #expected info is then block diagonal with one block for beta and one for tau:
  expected_info_half_inv <- Matrix::bdiag(expected_info_half_inv,
                                          Matrix::Diagonal(
                                            x = 1/sqrt(diag(expected_info)[tau_start:tau_end])))

  #compute augmentations
  augmentations <- matrix(0,nrow = n, ncol = J)

  for(i in 1:n){
    if(verbose){
      message("Computing data augmentation for observation ",i,".")
    }
    which_indices <- 1:J + (i - 1)*J #hmm maybe this isn't used...
    jacobian_i <- expected_info_list$jacobians[[i]] #jacobian (d[mean]/d(beta,tau), I believe)
    mu_i <- expected_info_list$mus[[i]] # means for i-th observation
    postmult_jacobian_i <- jacobian_i%*%  expected_info_half_inv #working toward asymmetric hat
    
    #entries of asymmetric hat matrix for observation i:
    H_i <- Matrix::tcrossprod(postmult_jacobian_i)%*%diag(as.numeric(mu_i))
    #diagonals are data augmentations:
    augmentations[i,] <- diag(as.matrix(H_i))[1:J]/2
  }

  return(augmentations)


}
