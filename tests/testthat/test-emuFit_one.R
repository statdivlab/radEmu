#
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
#   ml_fit <- emuFit_one(Y[,1],
#                        X = X,
#                        rect_weights = matrix(1,
#                                              ncol = J, nrow = n),
#                        j = 8
#                        B = rbind(b0,b1),
#                        z = z,
#                        method = "ML")
# })
