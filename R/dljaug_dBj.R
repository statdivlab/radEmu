
dljaug_dBj <- function(Y,
                       X,
                       B,
                       z,
                       constraint_fn,
                       constraint_grad_fn,
                       j,
                       j_constr,
                       k_constr
){
  mu_j <- exp(X%*%B[,j,drop = FALSE] + z)
  main_part <- t(X)%*%(Y[,j,drop = FALSE] - mu_j)

  mu_j_constr <- exp(X%*%B[,j_constr,drop = FALSE] + z)
  aug_scalar <- as.numeric(t(X[,k_constr,drop = FALSE])%*%
    (Y[,j_constr,drop = FALSE] - mu_j_constr))

  e_k_constr <- matrix(as.numeric((1:(ncol(X))) == k_constr),ncol = 1)

  dg <- constraint_grad_fn(B[k_constr,-j_constr])
  dg_index <- ifelse(j>j_constr,j - 1,j)
  dg <- dg[dg_index]
  aug_part <- e_k_constr*aug_scalar*dg # correct????

  return(main_part + aug_part)

}
