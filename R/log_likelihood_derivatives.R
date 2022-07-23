

log_likelihood_derivatives <- function(Y_tilde, X_tilde, log_mu_tilde,
                                       n, J,
                                       hessian = TRUE,
                                       A = NULL,
                                       Y_dot = NULL#,
                                       # loop_in_i = FALSE
                                       ){

#   if(loop_in_i){
#   i <- 1
#   ll_grad <- log_likelihood_derivatives(
#                         Y_tilde = Y_tilde[(i-1)*J + 1:J,,drop = FALSE],
#                         X_tilde = X_tilde[(i-1)*J + 1:J,],
#                         log_mu_tilde = log_mu_tilde[(i-1)*J + 1:J,,drop = FALSE],
#                         n = 1,
#                         J = J,
#                         A = A,
#                         hessian = FALSE,
#                         Y_dot = Y_dot[(i-1)*J + 1:J,,drop = FALSE],
#                         loop_in_i = FALSE)$gradient
#
#   for(i in 2:n){
#     ll_grad <- ll_grad +
#       log_likelihood_derivatives(
#       Y_tilde = Y_tilde[(i-1)*J + 1:J,,drop = FALSE],
#       X_tilde = X_tilde[(i-1)*J + 1:J,],
#       log_mu_tilde = log_mu_tilde[(i-1)*J + 1:J,,drop = FALSE],
#       n = 1,
#       J = J,
#       A = A,
#       hessian = FALSE,
#       Y_dot = Y_dot[(i-1)*J + 1:J,,drop = FALSE],
#       loop_in_i = FALSE)$gradient
#   }
#
# ll_hess <- NULL
# } else{

  jacob <- d_log_mu_tilde_d_B_tilde(Y_tilde, X_tilde, log_mu_tilde,n,J)

  if(is.null(A)){
  mini_A <- matrix(1,nrow = J, ncol = J)
  A <- Matrix::bdiag(lapply(1:n, function(i) mini_A))
  }
  if(is.null(Y_dot)){
  Y_dot <- A%*%Y_tilde
  }
  fitted_means <- Y_dot*exp(log_mu_tilde)

  resids <- Y_tilde - fitted_means

  ll_grad <- Matrix::crossprod(jacob,resids)


  # if(hessian){
  # ll_hess_pre <-
  #   Matrix::crossprod(jacob,
  #                     Matrix::Diagonal(n =length(fitted_means),
  #                                      x =sqrt(as.numeric(fitted_means))))
  #
  # ll_hess <- Matrix::tcrossprod(ll_hess_pre)
  # } else{
  #   ll_hess <- NULL
  # }
# }

  ll_hess <- NULL

  return(list("gradient" = ll_grad,
              "information" = ll_hess))
}
