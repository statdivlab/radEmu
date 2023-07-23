

constr_grad <- function(X,
                        Y,
                        B,
                        z,
                        j,
                        jstar,
                        kstar,
                        dg_dBj){



  if(j != jstar){
    mu_j = exp(X%*%B[,j,drop = FALSE] + z)
    mu_jstar = exp(X%*%B[,jstar,drop = FALSE] + z)
    G <- 0*X
    G[,kstar] = X[,kstar]*dg_dBj[kstar]

    dldBj = t(Y[,j,drop = FALSE] - mu_j)%*%X +
      t(Y[,jstar,drop = FALSE] - mu_jstar)%*%G
    return(dldBj)
  }

  if(j == jstar){
    mu_j = exp(X%*%B[,j,drop = FALSE] + z)
    dldBj = t(Y[,j,drop = FALSE] - mu_j)%*%X
    dldBj[,kstar] <- 0
    return(dldBj)
  }
}
