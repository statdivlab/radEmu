

calculate_log_likelihood <- function(Y_tilde, log_mu_tilde,
                                     J = NULL,
                                     n = NULL,
                                     A = NULL,
                                     Y_dot = NULL,
                                     log_Y_dot = NULL){

  if(is.null(A)){
    mini_A <- matrix(1,nrow = J, ncol = J)
    A <- Matrix::bdiag(lapply(1:n, function(i) mini_A))
  }
  if(is.null(Y_dot)){
  Y_dot <- A%*%Y_tilde
  }
  if(is.null(log_Y_dot)){
    log_Y_dot <- log(Y_dot)
  }

  return(sum(Y_tilde*(log_Y_dot + log_mu_tilde) - Y_dot*exp(log_mu_tilde)))
}
