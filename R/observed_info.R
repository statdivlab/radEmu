

observed_info <- function(X,Y,B){
  J <- ncol(Y)
  p <- ncol(X)
  n <- nrow(Y)
  z <- update_z_no_wts(Y,X,B)

  fitted_means <- exp(X%*%B + matrix(z,ncol = 1)%*%matrix(1,nrow = 1,ncol = J))
  observed_info <- matrix(0,ncol = p*J,nrow = p*J)
  X_iTX_is <- lapply(1:n,function(i) crossprod(X[i,,drop = FALSE]))
  p_denominators <- rowSums(exp(X%*%B))
  for(jprime in 1:J){
    for(jdoubleprime in jprime:J){
      jprime_indices <- (jprime - 1)*p + 1:p
      jdoubleprime_indices <- (jdoubleprime - 1)*p + 1:p

      p_jprimes <- exp(X%*%B[,jprime,drop = FALSE])/p_denominators

      p_jdoubleprimes <- exp(X%*%B[,jdoubleprime,drop = FALSE])/p_denominators

      for(j in 1:J){
        for(i in 1:n){
          observed_info[jprime_indices,jdoubleprime_indices] <-
            observed_info[jdoubleprime_indices,jprime_indices] <-
            observed_info[jdoubleprime_indices,jprime_indices] +
            (
              fitted_means[i,j]*
                (as.numeric(j == jdoubleprime) - p_jdoubleprimes[i])*
                  (as.numeric(j == jprime) - p_jprimes[i]) -
            (Y[i,j] - fitted_means[i,j])* p_jprimes[i]* p_jdoubleprimes[i])*
            X_iTX_is[[i]]
        }
      }
    }
  }

  return(observed_info)

}
