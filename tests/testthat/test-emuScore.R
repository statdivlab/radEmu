# test_that("Poisson unweighted CIs give reasonable results", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20))
#   z <- rnorm(40) +8
#   b0 <- rnorm(10)
#   b1 <- 0
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 10, nrow = 40)
#
#   for(i in 1:40){
#     for(j in 1:10){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rpois(1, lambda = temp_mean)
#     }
#   }
#   score_intervals <- emuScore(ml_fit)
#
#   score_intervals <- score_intervals[score_intervals$row == 2, ]
#
#   expect_true(mean(score_intervals$lower <= 0 & score_intervals$upper >=0) ==1)
# })
#
# test_that("Score intervals reasonable on NB data", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20))
#   z <- rnorm(40) +8
#   b0 <- rnorm(10)
#   b1 <- 0
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 10, nrow = 40)
#
#   for(i in 1:40){
#     for(j in 1:10){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rnbinom(1,mu = temp_mean,size = 1)
#     }
#   }
#   fl_fit <- emuFit_clone(Y = Y, X = X, method = "FL",
#                          constraint_fn = function(x) x[10])
#   score_intervals <- emuScore(fl_fit)
#
#   score_intervals <- score_intervals[score_intervals$row == 2, ]
#
#   expect_true(mean(score_intervals$lower <= 0 & score_intervals$upper >=0) ==1)
# })
#
# test_that("Score intervals from reweighted fit reasonable on NB data", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20))
#   z <- rnorm(40) +5
#   b0 <- rnorm(10)
#   b1 <- 0
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 10, nrow = 40)
#
#   for(i in 1:40){
#     for(j in 1:10){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rnbinom(1,mu = temp_mean,size = 1)
#     }
#   }
#   fl_fit <- emuFit_clone(Y = Y, X = X, method = "FL",
#                          reweight = TRUE,
#                          constraint_fn = function(x) x[10])
#   score_intervals <- emuScore(fl_fit)
#
#   score_intervals <- score_intervals[score_intervals$row == 2, ]
#
#   expect_equal(mean(score_intervals$lower <= 0 & score_intervals$upper >=0),
#               1,
#               tolerance = 0.1)
# })
