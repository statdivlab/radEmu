

profile_ll <- function(X,
                       Y,
                       B,
                       z = NULL){
  if(is.null(z)){
  z <- update_z_no_wts(Y,X,B)
  }
  log_mean <- X%*%B + matrix(z,ncol = 1)%*%matrix(1,nrow = 1, ncol = ncol(Y))
  return(sum(Y*log_mean - exp(log_mean)))
}
