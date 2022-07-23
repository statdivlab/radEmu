

X_to_X_tilde <- function(X,J){
  n <- nrow(X)
  p <- ncol(X)

  nonzero_i <- do.call(c,lapply(1:(n*J), function(i) rep(i,p)))
  nonzero_j <- do.call(c,lapply(1:n, function(i)
    do.call(c,lapply(1:J, function(k) 1:p + (k - 1)*p))))

  X_tilde <- Matrix::sparseMatrix(i = nonzero_i,
                                  j = nonzero_j,
                                  x = do.call(c,
                                              lapply(1:n,
                                                     function(i)
                                                       rep(X[i,],J))))
  return(X_tilde)

}
