
#convert B (square format) to B_cup (long format)
B_cup_from_B <- function(B){
  J <- ncol(B)
  p <- nrow(B)
  B_cup <- matrix(NA,nrow = J*p,ncol = 1)
  #so we just stack columns of B on top of each other to get B_cup
  for(j in 1:J){
    B_cup[(j - 1)*p + 1:p] <- B[,j]
  }
  return(Matrix::Matrix(B_cup))
}
