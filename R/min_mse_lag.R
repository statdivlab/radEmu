
min_mse_lag <- function(X,Y,B,constraint_grad_fn,k_constr,j_constr,j_ref){
  J <- ncol(Y)
  n <- nrow(Y)
  z <- update_z(Y,X,B)
  
  log_means <- X%*%B
  
  for(i in 1:n){
    log_means[i,] <- log_means[i,] + z[i]
  }
  

  
  d_dBjs <- lapply(1:J,
                  function(j)
                  Matrix::crossprod(X,
                                    Y[,j,drop = FALSE] - 
                                      exp(log_means[,j,drop = FALSE])))
  
  dldB <- do.call(cbind,d_dBjs)
  cg <- constraint_grad_fn(B[k_constr,-j_constr])[-(j_ref - as.numeric(j_ref > j_constr))]
  
  cg_mult <- dldB[k_constr,j_constr]
  
  cg <- cg*cg_mult
  
  lambda <- -sum(dldB[k_constr,-c(j_constr,j_ref)]*cg)/(sum(cg*cg))
  
  dldB[k_constr,-c(j_constr,j_ref)] <- dldB[k_constr,-c(j_constr,j_ref)]+
    lambda*cg

  dldB[k_constr,j_constr] <- 0
  dldB[,j_ref] <- 0
  
  return(sum(dldB^2))
}