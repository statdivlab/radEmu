

B_tilde_to_B <- function(B_tilde,J,p){

  B <- matrix(0,ncol = J,nrow = p)

  for(j in 1:J){
    B[,j] <-  B_tilde[(j - 1)*p + 1:p]
  }
  return(B)
}
