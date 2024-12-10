
#compute log likelihood after profiling out z
profile_ll <- function(X,Y,B){
  z <- update_z(Y = Y, X = X, B = B)
  log_means <- X%*%B + matrix(z,ncol = 1)%*%matrix(1,nrow = 1, ncol = ncol(Y))
  return(sum(Y*log_means - exp(log_means)))
}