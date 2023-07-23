
B_cup_from_B <- function(B){
  J <- ncol(B)
  p <- nrow(B)
  B_cup <- matrix(NA,nrow = J*p,ncol = 1)
  for(j in 1:J){
    B_cup[(j - 1)*p + 1:p] <- B[,j]
  }
  return(Matrix::Matrix(B_cup))
}
