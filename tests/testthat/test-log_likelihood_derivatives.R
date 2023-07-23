# test_that("Gradient of log likelihood matches numerical approximation", {
#
#
#
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20))
#   z <- rnorm(40) +8
#   b0 <- rnorm(10)
#   b1 <- 1:10
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 10, nrow = 40)
#
#   for(i in 1:40){
#     for(j in 1:10){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rpois(1, lambda = temp_mean)
#     }
#   }
#
#   n <- nrow(Y)
#   J <- ncol(Y)
#   p <- ncol(X)
#   B <- matrix(0,nrow = p, ncol = J)
#
#   get_ll <- function(B_tilde,J,n){
#     log_mu_tilde <- get_log_mu_tilde(X_tilde,B_tilde,n,J)
#     curr_ll <- calculate_log_likelihood(Y_tilde, log_mu_tilde, J, n)
#     return(curr_ll)
#   }
#
#   Y_tilde <- Y_to_Y_tilde(Y)
#   X_tilde <- X_to_X_tilde(X,J)
#   B_tilde <- B_to_B_tilde(B)
#
#   numerical_deriv <- numDeriv::grad(function(x)
#     get_ll(x,J = J, n = n),B_tilde)
#   log_mu_tilde <- get_log_mu_tilde(X_tilde,B_tilde, n, J)
#   analytical_deriv <- log_likelihood_derivatives(Y_tilde,
#                                                  X_tilde,
#                                                  log_mu_tilde,
#                                                  n,
#                                                  J,
#                                                  hessian = FALSE)
#
#   expect_equal(numerical_deriv,as.numeric(analytical_deriv$gradient))
#
# })
#
#
#
