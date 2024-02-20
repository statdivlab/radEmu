#get robust score stat
get_score_stat <- function(Y,
                           X_cup,
                           X,
                           B,
                           k_constr,
                           j_constr,
                           constraint_grad_fn,
                           indexes_to_remove, #indexes for B_vec corresp. to 
                                              #parameters set equal to zero as 
                                              #convenience constraint
                           j_ref, #index of taxon/category used as convenience constraint
                           J,
                           n,
                           p,
                           I_inv = NULL, #previously computed I_inv, if desired
                           Dy = NULL, #previously computed estimate of score covariance, if desired
                           cluster = NULL #numeric vector giving cluster membership for GEE
                           ){
  scores <- vector(n,mode = "list")
  
  p <- ncol(X)
  
  #reparametrize using convenience constraint
  for(k in 1:p){
    B[k,] <- B[k,] - B[k,j_ref]
  }
  
  #update z for new parametrization
  z <- update_z(Y = Y,
                X = X,
                B = B)
  
  #convert B to long format
  B_cup <- B_cup_from_B(B)
  #get score contribution of each observation
  for(i in 1:n) {
    
    X_cup_i <- X_cup[(i - 1)*J + 1:J, ]
    scores[[i]] <- Matrix::crossprod(X_cup_i, Y[i,] - exp(X_cup_i %*% B_cup + z[i]))
    
  }
  
  if(!is.null(cluster)){
    scores <- lapply(unique(cluster),
                     function(i) 
                       Reduce("+",scores[cluster == i]))
  }
  
  #compute derivative of constraint wrt (long/vector format) B 
  H <- matrix(0,nrow = p, ncol = J)
  
  H[k_constr,] <- constraint_grad_fn(B[k_constr,])
  
  H[k_constr,j_constr] <- H[k_constr,j_constr] - 1
  
  #we want derivative to be same dim as B_cup, so we convert:
  H_cup <- B_cup_from_B(H)
  
  #remove indexes corresponding to convenience constraint
  H_cup <- H_cup[-indexes_to_remove, , drop = FALSE]
  
  B_cup <- B_cup_from_B(B)
  
  if (is.null(I_inv)) {
    
    I <- f_info(Y, B_cup, B, X, X_cup)
    
    #remove indexes corresponding to convenience constraint
    I <- I[-indexes_to_remove, -indexes_to_remove]
    I_inv_H <- Matrix::solve(I, H_cup,  method = "cholmod_solve")
    
  } else {
    
    I_inv_H <- I_inv %*% H_cup
    
  }
  
  if (is.null(Dy)) {
    
    score_mat <- do.call(cbind, scores)
    score_mat <- methods::as(score_mat, "sparseMatrix")
    Dy <- Matrix::tcrossprod(score_mat)
    
  }
  #remove indexes corresponding to convenience constraint
  Dy <- Dy[-indexes_to_remove, -indexes_to_remove]
  
  
  score <- Reduce("+", scores)
  score <- as.matrix(score)
  #remove indexes corresponding to convenience constraint
  score <- score[-indexes_to_remove, , drop =FALSE]
  
  #adjustment factor from Guo et al. GEE paper (https://doi.org/10.1002/sim.2161)
  if(is.null(cluster)){
    score_adj <- n/(n - 1)
  } else{
    nclust <- length(unique(cluster))
    score_adj <- nclust/(nclust - 1)
  }
  
  #slightly fancy calculation of numerator and denominator of score stat:
  outside <- Matrix::crossprod(score, I_inv_H)
  inside <- Matrix::crossprod(I_inv_H, Dy) %*% I_inv_H
  score_stat <- as.numeric(as.matrix(outside^2/inside))*score_adj
  
  return(as.numeric(score_stat))
}
