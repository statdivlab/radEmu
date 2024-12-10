lil_ll <- function(shifts,B,p,X,Y,J,j_ref){
  temp_B <- B
  for(k in 1:p){
    temp_B[k,-j_ref] <-  temp_B[k,-j_ref] + shifts[k]
  }
  temp_z <- update_z(X = X, Y = Y, B = temp_B)
  log_means <- X%*%temp_B + matrix(temp_z, ncol = 1)%*%matrix(1,ncol = J,nrow = 1)
  
  return(-sum(Y*log_means - exp(log_means)))
}