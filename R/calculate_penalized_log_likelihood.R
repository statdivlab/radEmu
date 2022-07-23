# calculate_penalized_log_likelihood <- function(D_tilde_repar, theta_repar,
#                                      Y_tilde,
#                                      log_mu_tilde,
#                                      A = NULL,
#                                      Y_dot = NULL,
#                                      log_Y_dot = NULL){
#
#   if(is.null(A)){
#     mini_A <- matrix(1,nrow = J, ncol = J)
#     A <- Matrix::bdiag(lapply(1:n, function(i) mini_A))
#   }
#   if(is.null(Y_dot)){
#     Y_dot <- A%*%Y_tilde
#   }
#   if(is.null(log_Y_dot)){
#     log_Y_dot <- log(Y_dot)
#   }
#   info <- Matrix::crossprod(D_tilde,
#                             Matrix::crossprod(Diagonal(x = exp(as.numeric(log_mu_tilde))),D_tilde))
#   penalty <- 0.5*Matrix::determinant(info)$modulus
#
#   ll <- calculate_log_likelihood(Y_tilde,
#                                  log_mu_tilde,
#                                  A = A,
#                                  Y_dot = Y_dot,
#                                  log_Y_dot = log_Y_dot)
#
#   return(ll + penalty)
# }
