

calculate_firth_pl <-
  function(Y,
           B,
           X,
           D_tilde,
           z){
    J <- ncol(Y)
    n <- nrow(Y)
    B_tilde <- B_to_B_tilde(B)

    # X_tilde <- X_to_X_tilde(X,J = J)
    # dim(X_tilde)
    # nonzero_i <- rep(1:(n*J))
    # nonzero_j <- rep(1:n,each = J)
    # S_tilde <- Matrix::sparseMatrix(i = nonzero_i,
    #                                 j = nonzero_j,
    #                                 x = rep(1,n*J))
    # D_tilde <- cbind(X_tilde,S_tilde)
    W <- Matrix::Diagonal(x = exp(as.numeric(
      D_tilde%*%rbind(B_tilde,matrix(z,ncol = 1)))))
   penalty <- 0.5*Matrix::determinant(Matrix::crossprod(D_tilde,W)%*%D_tilde)$modulus

   logmeans <- X%*%B + matrix(z,ncol= 1)%*%matrix(1,nrow = 1, ncol = J)
   ll <- sum(Y*logmeans - exp(logmeans))

   return(as.numeric(ll + penalty))
  }
