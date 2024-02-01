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
                           Dy = NULL #previously computed estimate of score covariance, if desired
                           ){
  scores <- vector(n,mode = "list")
  
  p <- ncol(X)
  
  for(k in 1:p){
    B[k,] <- B[k,] - B[k,j_ref]
  }
  
  z <- update_z(Y = Y,
                X = X,
                B = B)
  
  B_cup <- B_cup_from_B(B)
  for(i in 1:n) {
    
    X_cup_i <- X_cup[(i - 1)*J + 1:J, ]
    scores[[i]] <- Matrix::crossprod(X_cup_i, Y[i,] - exp(X_cup_i %*% B_cup + z[i]))
    
  }
  
  H <- matrix(0,nrow = p, ncol = J)
  
  H[k_constr,] <- constraint_grad_fn(B[k_constr,])
  
  H[k_constr,j_constr] <- H[k_constr,j_constr] - 1
  
  
  H_cup <- B_cup_from_B(H)
  
  H_cup <- H_cup[-indexes_to_remove, , drop = FALSE]
  
  B_cup <- B_cup_from_B(B)
  
  if (is.null(I_inv)) {
    
    I <- f_info(Y, B_cup, B, X, X_cup)
    
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
  
  Dy <- Dy[-indexes_to_remove, -indexes_to_remove]
  
  
  score <- Reduce("+", scores)
  score <- as.matrix(score)
  score <- score[-indexes_to_remove, , drop =FALSE]
  
  #slightly fancy calculation of numerator and denominator of score stat:
  outside <- Matrix::crossprod(score, I_inv_H)
  inside <- Matrix::crossprod(I_inv_H, Dy) %*% I_inv_H
  score_stat <- as.numeric(as.matrix(outside^2/inside))
  
  return(as.numeric(score_stat))
}
