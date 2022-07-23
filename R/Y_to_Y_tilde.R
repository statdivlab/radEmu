

Y_to_Y_tilde <- function(Y){

  n <- nrow(Y)
  J <- ncol(Y)
  Y_tilde <- matrix(0, ncol =1 , nrow = n*J)

  for(i in 1:n){
    Y_tilde[(i - 1)*J + 1:J] <- Y[i,]
  }
  return(Y_tilde)
}
