

get_C <- function(A,
                  null_index){
  A21 <- A[null_index,-null_index]
  A11 <- A[-null_index,-null_index]

  A11_ginv <- MASS::ginv(as.matrix(A11))


  C <- cbind(-A21%*%A11_ginv,1)
  return(C)

}
