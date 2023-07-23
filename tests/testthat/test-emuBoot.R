# test_that("emuBoot returns appropriately formatted output", {
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
#   ml_fit <- fitted_model_ml<- emuFit(X = X, Y = Y,method = "ML",
#                                      tolerance = 0.01)
#
#   boot_runs <- emuBoot(ml_fit,
#                        10)
#   expect_true(length(boot_runs)==10)
#   expect_true(unique(sapply(1:10,function(x) nrow(boot_runs[[x]]))) ==2)
#   expect_true(unique(sapply(1:10,function(x) ncol(boot_runs[[x]]))) ==10)
# })
