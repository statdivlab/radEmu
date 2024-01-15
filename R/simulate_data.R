#' @importFrom stats rnorm rbinom rpois rnbinom

simulate_data <- function(n,
                          J,
                          b0,
                          b1,
                          distn,
                          zinb_size = NULL,
                          zinb_zero_prop = NULL,
                          mean_count_before_ZI) {
  
  X <- cbind(1,rep(c(0,1),each = n/2))
  B <- rbind(b0,b1)
  log_means <- do.call(cbind,
                       lapply(1:J,
                              function(j) X%*%B[,j,drop = FALSE]))

  row_means <- rowSums(exp(log_means))/J

  z <- sapply(row_means,function(x) log(mean_count_before_ZI) - log(x) + stats::rnorm(1))
  Y <- matrix(0, ncol = J, nrow = n)

  for(i in 1:n){
    log_means[i,] <- log_means[i,] + z[i]
  }

  for(i in 1:n){
    accepted <- FALSE
    while(!accepted){
    for(j in 1:J){
      if(distn == "Poisson"){
        Y[i,j] <- stats::rpois(1,lambda = exp(log_means[i,j]))
      }
      if(distn == "ZINB"){
        Y[i,j] <- stats::rnbinom(1,mu = exp(log_means[i,j]), size= zinb_size)*
          (1 - stats::rbinom(1,1,prob = zinb_zero_prop))
      }
      if(sum(Y[i,])>0){
        accepted <- TRUE
      }
      }
    }
  }
  return(Y)
}
