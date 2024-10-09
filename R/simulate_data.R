#' Data simulation function
#' 
#' Function to simulate data for simulations in Clausen & Willis (2024) and for the cluster vignette
#'
#' @param n Number of samples
#' @param J Number of categories 
#' @param b0 Intercept parameter vector
#' @param b1 Covariate paramter vector 
#' @param distn Distribution to simulate from, either "Poisson" or "ZINB"
#' @param zinb_size Size parameter for negative binomial draw for ZINB data
#' @param zinb_zero_prop Proportion of zeros for ZINB data
#' @param mean_count_before_ZI Parameter for construction of z_i in mean model. Setting this to \code{50} works well in practice. 
#' @param X Optional design matrix, this must have two columns and n rows.
#' @param B Optional B matrix, if p is not equal to 2 
#' @param cluster Optional cluster vector, this must have n elements. 
#'
#' @return \code{Y}. A \code{n times J} dimension matrix of simulated response counts.
#'
#' @importFrom stats rnorm rbinom rpois rnbinom
#'
#' @export
simulate_data <- function(n,
                          J,
                          b0 = NULL,
                          b1 = NULL,
                          distn,
                          zinb_size = NULL,
                          zinb_zero_prop = NULL,
                          mean_count_before_ZI,
                          X = NULL,
                          B = NULL,
                          cluster = NULL) {
  
  if (is.null(X)) {
    X <- cbind(1,rep(c(0,1),each = n/2))
  }
  if (!is.null(b0) & !is.null(b1)) {
    B <- rbind(b0,b1)
    if (nrow(B) != ncol(X)) {
      stop("You've input b0 and b1 but your X matrix does not have 2 columns. Please use the B argument when your design matrix does not have 2 columns.")
    }
  } else if (is.null(b0) | is.null(b1)) {
    if (is.null(B)) {
      stop("Please input either parameter vectors b0 and b1, or parameter matrix B.")
    }
  }
  
  log_means <- do.call(cbind,
                       lapply(1:J,
                              function(j) X %*% B[,j,drop = FALSE]))
  
  z <- mean_count_before_ZI + stats::rnorm(n)
  
  Y <- matrix(0, ncol = J, nrow = n)

  for(i in 1:n){
    log_means[i,] <- log_means[i,] + z[i]
  }
  
  if (is.null(cluster)) {
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
  } else {
    cluster_effs <- lapply(1:length(unique(cluster)), function(i) log(matrix(stats::rexp(2*J), nrow = 2)))
    for(i in 1:n){
      accepted <- FALSE
      while(!accepted){
        for(j in 1:J){
          if(distn == "Poisson"){
            Y[i,j] <- stats::rpois(1,lambda = exp(log_means[i,j] + 
                                                    cluster_effs[[cluster[i]]][, j]))
          }
          if(distn == "ZINB"){
            Y[i,j] <- stats::rnbinom(1, 
                                     mu = exp(log_means[i,j] + 
                                                cluster_effs[[cluster[i]]][, j]), 
                                     size= zinb_size)*
              (1 - stats::rbinom(1,1,prob = zinb_zero_prop))
          }
          if(sum(Y[i,])>0){
            accepted <- TRUE
          }
        }
      }
    }
  }
  
  return(Y)
}
