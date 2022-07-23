

log_likelihood_wide <- function(Y,
                                rect_weights,
                                X,
                                B,
                                z){
  J <- ncol(Y)
  log_means <- X%*%B + matrix(z,ncol = 1)%*%matrix(1,ncol = J,nrow = 1)
  return(sum((Y*log_means - exp(log_means))*rect_weights))
}
