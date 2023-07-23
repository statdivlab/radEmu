Y_long_from_Y <- function(Y){
  n <- nrow(Y)
  J <- ncol(Y)
  Y_long <- matrix(nrow = n*J,ncol = 1)

  for(i in 1:n){
    Y_long[1:J + J*(i - 1),] <- Y[i,]
  }

  return(Y_long)
}
