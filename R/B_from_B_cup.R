

B_from_B_cup <- function(B_cup,J,p){
  B <- matrix(0,nrow = p, ncol = J)

  for(j in 1:J){
    B[,j] <-  B_cup[(j - 1)*p + 1:p]
  }
  return(B)
}
