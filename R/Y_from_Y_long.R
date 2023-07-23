Y_from_Y_long <- function(Y_long, n, J){
  Y <- matrix(nrow = n, ncol = J)
  for(i in 1:n){
    Y[i,] <- Y_long[1:J + (i - 1)*J,]
  }
  return(Y)
}
