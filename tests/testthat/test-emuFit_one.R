# #
# test_that("ML fit to simple example give reasonable output", {
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
#   ml_fit <- emuFit_one(Y = Y,
#                        X = X,
#                        rect_weights = 0*Y + 1,
#                        j = 8,
#                        B = rbind(b0,b1),
#                        z = z,
#                        method = "ML",
#                        maxit_glm = 50,
#                        stop_on_error = FALSE,
#                        WD = NULL,
#                        info_inv = NULL,
#                        constraint_fn = NULL,
#                        constrained_j = NULL,
#                        constrained_k = NULL)
#
#   expect_true(inherits(ml_fit,"numeric"))
#   B[,8] <- ml_fit
#
#
# })
#
#
# #
# test_that("ML fit to simple example give reasonable output with constraints", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20))
#   z <- rnorm(40) +8
#   b0 <- rnorm(10)
#   b0 <- b0 - pseudohuber_center(b0)
#   b1 <- 1:10
#   b1 <- b1 - pseudohuber_center(b1[-2])
#   b1[2] <- 0
#   b <- rbind(b0,b1)
#
#   Y <- matrix(NA,ncol = 10, nrow = 40)
#
#   for(i in 1:40){
#     for(j in 1:10){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rpois(1, lambda = temp_mean)
#     }
#   }
#
#   ml_fit <- emuFit_one(Y = Y,
#                        X = X,
#                        rect_weights = 0*Y + 1,
#                        j = 8,
#                        B = b,
#                        z = z,
#                        method = "ML",
#                        maxit_glm = 50,
#                        stop_on_error = FALSE,
#                        WD = NULL,
#                        info_inv = NULL,
#                        constraint_fn = pseudohuber_center,
#                        constrained_j = 2,
#                        constrained_k = 2)
#
#   expect_true(inherits(ml_fit,"numeric"))
#   b[,8] <- ml_fit
#
#   expect_equal(pseudohuber_center(b[2,]),as.numeric(b[2,2]),
#                tolerance = 1e-4)
#
#
# })
#
