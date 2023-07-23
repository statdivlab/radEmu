

constr_info <- function(X,
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
    I <- t(X)%*%diag(as.numeric(mu_j))%*%X  +
      t(G)%*%diag(as.numeric(mu_jstar))%*%G
    return(I)
  }
  if(j == jstar){
    mu_j = exp(X%*%B[,j,drop = FALSE] + z)
    Xstar <- X
    Xstar[,kstar] <- 0*Xstar[,kstar]
    I <- t(Xstar)%*%diag(as.numeric(mu_j))%*%Xstar
    return(I)
  }

  }
