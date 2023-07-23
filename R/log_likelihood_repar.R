

log_likelihood_repar <- function(Y_long,
                                 beta_tilde_J,
                                 tau,
                                 X_tilde_J,
                                 J,
                                 n,
                                 p){

  log_means <- X_tilde_J%*%beta_tilde_J
  for(i in 1:n){
    log_means[1:J + (i - 1)*J,,drop = FALSE] <-
      log_means[1:J + (i - 1)*J,,drop = FALSE] -
      log(sum(exp(log_means[1:J + (i - 1)*J,,drop = FALSE]))) + log(tau[i])
  }

  ll <- sum(Y_long*log_means - exp(log_means))

  return(ll)
}
