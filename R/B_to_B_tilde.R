

B_to_B_tilde <- function(B){

  p <- nrow(B)
  J <- ncol(B)
  B_tilde <- matrix(0, ncol =1 , nrow = p*J)
  #
  # for(k in 1:p){
  #   B_tilde[(k - 1)*J + 1:J] <- B[k,]
  # }
  # return(B_tilde)
  for(j in 1:J){
    B_tilde[(j - 1)*p + 1:p] <- B[,j]
  }
  return(B_tilde)
}
