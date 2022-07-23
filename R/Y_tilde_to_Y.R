
Y_tilde_to_Y <- function(Y_tilde,J){
  if(!is.matrix(Y_tilde)){
    Y_tilde <- matrix(as.numeric(Y_tilde),ncol = 1)
  }
  n <- nrow(Y_tilde)/J
  if((n - round(n,0)) != 0){
    stop("Number of entries in Y_tilde must be a multiple of J.")
  }
  Y <- matrix(0,ncol = J, n)

  for(i in 1:n){
    Y[i,] <-  Y_tilde[(i - 1)*J + 1:J]
  }
return(Y)
}
