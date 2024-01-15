# test_that("Analytical hessian of Lagrangian equal to numerical hessian", {
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
#   u <- 1e8
#   rho <- 1e7
#
#
#   num_constr_hess <- numDeriv::hessian(function(x) pseudohuber_center(x,0.1),
#                                        x = ml_fit[2,],
#                                        method = "Richardson",
#                                        method.args=list(eps=1e-2))
#
#   an_constr_hess <- matrix(NA, nrow = 10, ncol = 10)
#
#   for(i in 1:10){
#     for(j in i:10){
#       an_constr_hess[i,j] <- an_constr_hess[j,i] <-
#         hess_pseudohuber_center(ml_fit[2,],0.1,i,j)
#     }
#   }
#
#   # an_constr_hess/ num_constr_hess
#
#   another_num_hess <- do.call(cbind,lapply(1:10, function(j)
#                                            numDeriv::grad(function(x) dpseudohuber_center_dx(x,0.1)[j],
#                                      x = ml_fit[2,])))
#
#   # max(abs(an_constr_hess/another_num_hess - 1))
#
#   expect_equal(an_constr_hess,
#                another_num_hess,tolerance = 1e-3)
#
#
#   # message("test that numerical derivatives of lagrangian_grad agree with analytical hessian")
#   an_hess <- lagrangian_hessian(X = X,
#                                 Y = Y,
#                                 B = ml_fit,
#                                 z = z,
#                                 j = 1,
#                                 j_constr = 3,
#                                 k_constr = 2,
#                                 u = u,
#                                 rho = rho,
#                                 feas_gap = ml_fit[2,3] - pseudohuber_center(ml_fit[2,]),
#                                 constraint_grad_fn = function(x) dpseudohuber_center_dx(x,0.1),
#                                 constraint_hess_fn = function(x,ind_1,ind_2)
#                                   hess_pseudohuber_center(x,0.1,ind_1,ind_2))
#
#   num_hess <- numDeriv::hessian(function(x){
#     temp_B <- ml_fit
#     temp_B[,1] <- x
#     feas_gap <- temp_B[2,3] - pseudohuber_center(temp_B[2,],0.1)
#     return(  lagrangian(X,
#                         Y,
#                         B = temp_B,
#                         z,
#                         u = u,
#                         k_constr = 2,
#                         j_constr = 3,
#                         rho = rho,
#                         feas_gap = feas_gap))
#   },
#   ml_fit[,1])
#
#   another_num_hess <-
#   do.call(cbind,
#           lapply(1:2,
#                  function(kay)
#                    numDeriv::grad(func = function(bj){
#                      temp_B <- ml_fit
#                      temp_B[,1] <- bj
#                      feas_gap <- temp_B[2,3] - pseudohuber_center(temp_B[2,],0.1)
#                      lagrangian_grad(X,
#                                      Y,
#                                      B = temp_B,
#                                      z,
#                                      u = u,
#                                      k_constr = 2,
#                                      j_constr = 3,
#                                      rho = rho,
#                                      constraint_grad_fn =
#                                        function(w) dpseudohuber_center_dx(w,0.1),
#                                      feas_gap = feas_gap)[[1]][kay]
#                    },
#                      x = ml_fit[,1])
#           ))
#
#
#
#
#
#   # an_hess/num_hess
#   # num_hess/another_num_hess
#   expect_equal(num_hess,an_hess,tolerance = 1e-3)
#
# })
#
#
# test_that("Lagrangian hessian function warns you if constraint_hess_fn is not provided", {
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
#   expect_warning(lagrangian_hessian(X = X,
#                                 Y = Y,
#                                 B = ml_fit,
#                                 z = z,
#                                 j = 1,
#                                 j_constr = 3,
#                                 k_constr = 2,
#                                 u = 1,
#                                 rho = 1,
#                                 feas_gap = ml_fit[2,3] - pseudohuber_center(ml_fit[2,]),
#                                 constraint_grad_fn = function(x) dpseudohuber_center_dx(x,0.1)))
#
# })
#
