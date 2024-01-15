


update_z <- function(Y,
                     X,
                     B){

    z <- log(Matrix::rowSums(Y)) -
      log(Matrix::rowSums(exp(X%*%B)))

  return(z)
}
