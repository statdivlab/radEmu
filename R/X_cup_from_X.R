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

X_cup_from_X_fast <- function(X, J) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Total number of rows in final matrix
  total_rows <- n * J
  
  # Construct i index (row indices)
  i_coords <- rep(1:total_rows, each = p)
  
  # Construct j index (column indices)
  j_block <- rep(0:(J - 1), each = p) * p
  j_coords <- rep(j_block, times = n) + rep(1:p, times = total_rows)
  
  # Construct values
  X_rep <- X[rep(1:n, each = J), ]
  x_vals <- as.vector(t(X_rep))  # row-wise unrolling
  
  X_cup <- Matrix::sparseMatrix(i = i_coords,
                                j = j_coords,
                                x = x_vals)
  
  return(X_cup)
}
