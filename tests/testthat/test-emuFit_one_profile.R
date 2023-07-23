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
#
#  full_fit <-  emuFit_micro(Y = Y,
#          X = X,
#          constraint_fn = pseudohuber_center
#          )
#
#   for(iter in 1:100){
#     print(iter)
#     for(j in 1:10){
#   ml_fit <- emuFit_one_profile(Y = Y,
#                        X = X,
#                        j = 8,
#                        B - b,
#                        maxit_glm = 50,
#                        constraint_fn = NULL,
#                        constrained_j = NULL,
#                        constrained_k = NULL)
#
#   b[,j] <- ml_fit
#
#     }
#   b[1,] <- b[1,] - median(b[1,])
#   b[2,] <- b[2,] - median(b[2,])
#     print(b)
#   }
#
#   plot(b[1,])
#   plot(b[2,])
#
#   expect_true(inherits(ml_fit,"numeric"))
#   B[,8] <- ml_fit
#
#
# })
