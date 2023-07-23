X_cup_from_X <- function(X,J){
  p <- ncol(X)

  n <- nrow(X)

  i_coords <- numeric(0)
  j_coords <- numeric(0)

  for(i in 1:n){
    i_coords <- c(i_coords,
                  rep((i - 1)*J + rep(1:J, each = p)))
    j_coords <- c(j_coords,
                  sapply(1:J, function(j) (j - 1)*p + 1:p))
  }

  X_cup <- Matrix::sparseMatrix(x = do.call(c,lapply(1:n,
                                            function(i) rep(X[i,],J))),
                                i = i_coords,
                                j = j_coords)

  return(X_cup)
}
