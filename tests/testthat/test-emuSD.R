# test_that("emuSD returns appropriately formatted output", {
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
#   boot_sd <- emuSD(ml_fit,
#                    nboot = 10,
#                    rows_of_interest = c(1,2))
#   expect_true(nrow(boot_sd) ==2)
#   expect_true(ncol(boot_sd) ==10)
#   expect_true(sum(boot_sd <=0) ==0)
# })
