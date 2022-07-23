


update_z <- function(Y,
                     rect_weights,
                     X,
                     B){
    z <- log(Matrix::rowSums(Y*rect_weights)) -
      log(Matrix::rowSums(exp(X%*%B)*rect_weights))
    return(z)
}
