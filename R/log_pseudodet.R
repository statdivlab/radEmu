log_pseudodet <- function(A,mat_rank,delta = 1e-5){

  if(ncol(A)!=nrow(A)){
    stop("Matrix A must be square")
  }

  mat_n <- nrow(A)
  mat_k <- mat_rank


  A_prime <- A + diag(rep(delta,ncol(A)))

  return(Matrix::determinant(A_prime)$modulus - (mat_n - mat_k)*log(delta))
}
