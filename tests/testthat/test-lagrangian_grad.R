# test_that("multiplication works", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20))
#   z <- rnorm(40) +5
#   J <- 10
#   p <- 2
#   n <- 40
#   b0 <- rnorm(J)
#   b1 <- seq(1,10,length.out = J)
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = J, nrow = 40)
#
#   for(i in 1:40){
#     for(j in 1:J){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rpois(1, lambda = temp_mean)
#     }
#   }
#   ml_fit <- emuFit_micro(X,
#                          Y,
#                          maxit = 200,
#                          tolerance = 1e-2)
#
#   z <- update_z(Y,X,ml_fit)
#   lag_grad <-
#     lagrangian_grad(X = X,
#                   Y = Y,
#                   B = ml_fit,
#                   z = z,
#                   u = 1e6,
#                   rho = 1e9,
#                   j = 1,
#                   k_constr = 2,
#                   j_constr = 3,
#                   feas_gap = ml_fit[2,2] - pseudohuber_center(ml_fit[2,],0.1),
#                   constraint_grad_fn = function(x){
#                     dpseudohuber_center_dx(x,0.1)
#                   })
#
#   num_grad <- numDeriv::grad(function(x){
#     temp_B <- ml_fit
#     temp_B[,1] <- x
#     return(
#       lagrangian(X = X,
#                  Y = Y,
#                  B = temp_B,
#                  z = z,
#                  u = 1e6,
#                  rho = 1e9,
#                  k_constr = 2,
#                  j_constr = 3,
#                  feas_gap = temp_B[2,2] - pseudohuber_center(temp_B[2,],0.1)
#                  )
#     )
#   },
#   x = ml_fit[,1])
#
#   expect_equal(lag_grad[[1]],num_grad,tolerance = 1e-3)
#
#   lag_grad <-
#     lagrangian_grad(X = X,
#                     Y = Y,
#                     B = ml_fit,
#                     z = z,
#                     u = -1e7,
#                     rho = 1e6,
#                     j = 1,
#                     k_constr = 2,
#                     j_constr = 8,
#                     feas_gap = ml_fit[2,8] - pseudohuber_center(ml_fit[2,],0.1),
#                     constraint_grad_fn = function(x){
#                       dpseudohuber_center_dx(x,0.1)
#                     })
#
#   num_grad <- numDeriv::grad(function(x){
#     temp_B <- ml_fit
#     temp_B[,1] <- x
#     return(
#       lagrangian(X = X,
#                  Y = Y,
#                  B = temp_B,
#                  z = z,
#                  u = -1e7,
#                  rho = 1e6,
#                  k_constr = 2,
#                  j_constr = 8,
#                  feas_gap = temp_B[2,8] - pseudohuber_center(temp_B[2,],0.1)
#       )
#     )
#   },
#   x = ml_fit[,1])
#
#   expect_equal(lag_grad[[1]],num_grad,tolerance = 1e-3)
# })
