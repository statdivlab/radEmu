# test_that("emuFit_clone works with weights = 1", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20))
#   z <- rnorm(40) + 8
#   b0 <- rnorm(10)
#   b1 <- 1:10
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 10, nrow = 40)
#
#   for(i in 1:40){
#     for(j in 1:10){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rnbinom(1, mu = temp_mean,size = 1)
#     }
#   }
#   weights <- matrix(rep(1,40),ncol = 1)%*%matrix(1,ncol = 10,
#                                                  nrow = 1)
#
#   non_cheat_fit <- emuFit(X = X, Y = Y,method = "ML",
#                           fast_fit = FALSE,
#                           fast_fit_cheat = FALSE,
#                           verbose = TRUE,
#                           reweight = FALSE,
#                           weights = weights,
#                           tolerance = 0.01)
#
#   clone_fit <- emuFit_clone(X = X, Y = Y,
#                             method = "ML",
#                             weights = weights)
#
#   expect_equal(clone_fit$B,non_cheat_fit$B,tolerance = 1e-6)
# })
#
# test_that("emuFit_clone works when n >>j", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 200))
#   z <- rnorm(400) +8
#   b0 <- rnorm(10)
#   b1 <- 1:10
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 10, nrow = 400)
#
#   for(i in 1:400){
#     for(j in 1:10){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rnbinom(1, mu = temp_mean,size = 0.25)
#     }
#   }
#   weights <- matrix(rep(1,400),ncol = 1)%*%matrix(1,ncol = 10,
#                                                   nrow = 1)
#
#   non_cheat_fit <- emuFit(X = X, Y = Y,method = "ML",
#                           fast_fit = FALSE,
#                           fast_fit_cheat = FALSE,
#                           verbose = TRUE,
#                           reweight = FALSE,
#                           weights = weights,
#                           tolerance = 0.01)
#
#   clone_fit <- emuFit_clone(X = X, Y = Y,
#                             method = "ML",
#                             weights = weights,
#                             tolerance = 0.01)
#
#   expect_equal(clone_fit$B,non_cheat_fit$B,tolerance = 1e-3)
# })
#
# test_that("emuFit_clone works with Firth likelihood", {
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
#       Y[i,j] <- rnbinom(1, mu = temp_mean,size = 0.25)
#     }
#   }
#   weights <- matrix(rep(1,40),ncol = 1)%*%matrix(1,ncol = 10,
#                                                  nrow = 1)
#
#   non_cheat_fit <- emuFit(X = X, Y = Y,method = "FL",
#                           fast_fit = FALSE,
#                           fast_fit_cheat = FALSE,
#                           verbose = TRUE,
#                           reweight = FALSE,
#                           weights = weights,
#                           tolerance = 0.01)
#
#   clone_fit <- emuFit_clone(X = X, Y = Y,
#                             method = "FL",
#                             weights = weights)
#
#   expect_equal(clone_fit$B,non_cheat_fit$B,tolerance = 1e-5)
# })
#
#
# test_that("emuFit_clone works when n << j", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20))
#   z <- rnorm(40) +8
#   b0 <- rnorm(200)
#   b1 <- (1:200)/20
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 200, nrow = 40)
#
#   for(i in 1:40){
#     for(j in 1:200){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rnbinom(1, mu = temp_mean,size = 0.25)
#     }
#   }
#   weights <- matrix(rep(1,40),ncol = 1)%*%matrix(1,ncol = 200,
#                                                  nrow = 1)
#
#   non_cheat_fit <- emuFit(X = X, Y = Y,method = "ML",
#                           fast_fit = FALSE,
#                           fast_fit_cheat = FALSE,
#                           verbose = TRUE,
#                           reweight = FALSE,
#                           weights = weights,
#                           tolerance = 0.001)
#
#   clone_fit <- emuFit_clone(X = X, Y = Y,
#                             method = "ML",
#                             weights = weights)
#
#   expect_equal(clone_fit$B,non_cheat_fit$B,tolerance = 1e-3)
# })
#
# test_that("emuFit_clone works with Firth likelihood when n << j", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20))
#   z <- rnorm(40) +8
#   b0 <- rnorm(200)
#   b1 <- (1:200)/20
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 200, nrow = 40)
#
#   for(i in 1:40){
#     for(j in 1:200){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rnbinom(1, mu = temp_mean,size = 0.25)
#     }
#   }
#   weights <- matrix(rep(1,40),ncol = 1)%*%matrix(1,ncol = 200,
#                                                  nrow = 1)
#
#   non_cheat_fit <- emuFit(X = X, Y = Y,method = "FL",
#                           fast_fit = FALSE,
#                           fast_fit_cheat = FALSE,
#                           verbose = TRUE,
#                           reweight = FALSE,
#                           weights = weights,
#                           tolerance = 0.001)
#
#   clone_fit <- emuFit_clone(X = X, Y = Y,
#                             method = "FL",
#                             weights = weights)
#
#   expect_equal(clone_fit$B,non_cheat_fit$B,tolerance = 1e-5)
# })
#
# test_that("emuFit_clone works with nb data", {
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
#       Y[i,j] <- rnbinom(1, mu = temp_mean,size = 0.1)
#     }
#   }
#   weights <- matrix(rep(1,40),ncol = 1)%*%matrix(1,ncol = 10,
#                                                  nrow = 1)
#
#   non_cheat_fit <- emuFit(X = X, Y = Y,method = "ML",
#                           fast_fit = FALSE,
#                           fast_fit_cheat = FALSE,
#                           verbose = TRUE,
#                           reweight = FALSE,
#                           weights = weights,
#                           tolerance = 0.01)
#
#   clone_fit <- emuFit_clone(X = X, Y = Y,
#                             method = "ML",
#                             weights = weights,
#                             tolerance = 0.01)
#
#   expect_equal(clone_fit$B,non_cheat_fit$B,tolerance = 1e-3)
# })
#
#
# test_that("emuFit_clone works with weights", {
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
#       Y[i,j] <- rnbinom(1, mu = temp_mean,size = 0.25)
#     }
#   }
#   weights <- matrix(rgamma(40,shape= 1/10,rate = 1),
#                     ncol = 1)%*%matrix(1,ncol = 10,
#                                        nrow = 1)
#
#   non_cheat_fit <- emuFit(X = X, Y = Y,method = "ML",
#                           fast_fit = FALSE,
#                           fast_fit_cheat = FALSE,
#                           verbose = TRUE,
#                           reweight = FALSE,
#                           weights = weights,
#                           tolerance = 0.01)
#
#   clone_fit <- emuFit_clone(X = X, Y = Y,
#                             method = "ML",
#                             weights = weights)
#
#   expect_equal(clone_fit$B,non_cheat_fit$B,tolerance = 1e-5)
# })
#
#
# test_that("emuFit_clone works with weights differing in j", {
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
#       Y[i,j] <- rnbinom(1, mu = temp_mean,size = 0.25)
#     }
#   }
#   weights <- matrix(1,ncol = 1, nrow = 40)%*%matrix(1:10,
#                                                     ncol = 10, nrow = 1)
#
#   non_cheat_fit <- emuFit(X = X, Y = Y,method = "ML",
#                           fast_fit = FALSE,
#                           fast_fit_cheat = FALSE,
#                           verbose = TRUE,
#                           reweight = FALSE,
#                           weights = weights,
#                           tolerance = 0.01)
#
#   clone_fit <- emuFit_clone(X = X, Y = Y,
#                             method = "ML",
#                             weights = weights)
#
#   expect_equal(clone_fit$B,non_cheat_fit$B,tolerance = 1e-4)
# })
#
#
# test_that("emuFit_clone works with weights correlated with X", {
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
#       # Y[i,j] <- rnbinom(1, mu = temp_mean,size = 0.25)
#       Y[i,j] <- rpois(1,temp_mean)
#     }
#   }
#   weights <- matrix(c(rep(1,20),rep(10,20)),ncol = 1)%*%matrix(1,ncol =10, nrow = 1)
#
#   non_cheat_fit <- emuFit(X = X, Y = Y,method = "ML",
#                           fast_fit = FALSE,
#                           fast_fit_cheat = FALSE,
#                           verbose = TRUE,
#                           reweight = FALSE,
#                           weights = weights,
#                           tolerance = 0.0001)
#
#   clone_fit <- emuFit_clone(X = X, Y = Y,
#                             method = "ML",
#                             weights = weights)
#
#   expect_equal(clone_fit$B,non_cheat_fit$B,tolerance = 1e-1)
# })
#
# test_that("emuFit_clone with weights returns same result as orig fit", {
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
#   weights <- matrix(rgamma(40,shape= 1/100,rate = 1),
#                     ncol = 1)%*%matrix(1,ncol = 10,
#                                        nrow = 1)
#
#   clone_fit <- emuFit_clone(Y = Y, X = X,
#                             method = "ML",
#                             # reweight = TRUE
#                             weights = weights)
#
#
#   clone_fitted <-
#     X%*%clone_fit$B +
#     matrix(clone_fit$z,ncol = 1)%*%matrix(1,nrow = 1,ncol = 10)
#
#
#   non_cheat_fit <- emuFit(X = X, Y = Y,method = "ML",
#                           fast_fit = FALSE,
#                           fast_fit_cheat = FALSE,
#                           verbose = TRUE,
#                           # reweight = TRUE,
#                           weights= weights,
#                           tolerance = 1e-10)
#
#   non_cheat_fitted <-
#     X%*%non_cheat_fit$B +
#     matrix(non_cheat_fit$z,ncol = 1)%*%matrix(1,nrow = 1,ncol = 10)
#
#   expect_equal(non_cheat_fit$ll, clone_fit$ll[[2]])
#   expect_true(max(abs(clone_fit$B[2,] - seq(-4.5,4.5,1)))<.1)
#   expect_equal(clone_fit$B,non_cheat_fit$B,tolerance = 1e-5)
#   expect_equal(clone_fit$weights,non_cheat_fit$weights)
# })
#
# test_that("emuFit_clone with reweighting returns same result as orig fit", {
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
#   clone_fit <- emuFit_clone(Y = Y, X = X,
#                             method = "ML",
#                             reweight = TRUE)
#
#
#
#   non_cheat_fit <- emuFit(X = X, Y = Y,method = "ML",
#                           fast_fit = FALSE,
#                           fast_fit_cheat = FALSE,
#                           verbose = TRUE,
#                           reweight = TRUE,
#                           # weights= weights,
#                           tolerance = 1e-20)
#
#
#   clone_fit_preweight <-
#     emuFit_clone(Y = Y, X = X,
#                  method = "ML",
#                  weights = clone_fit$weights,
#                  verbose = TRUE)
#
#
#   clone_preweight_fitted <-
#     X%*%clone_fit_preweight$B +
#     matrix(clone_fit_preweight$z,ncol = 1)%*%matrix(1,nrow = 1,ncol = 10)
#
#   expect_equal(clone_fit$B, clone_fit_preweight$B,tolerance = 1e-5)
#   expect_equal(non_cheat_fit$B, clone_fit_preweight$B,tolerance = 1e-1)
#   expect_equal(clone_fit$B,clone_fit_preweight$B,tolerance = 1e-4)
#   expect_equal(non_cheat_fit$B, clone_fit$B,tolerance= 1e-1)
#
# })
#
# test_that("ML emuFit_clone returns same result as orig fit on sparse Y", {
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
#       Y[i,j] <- rpois(1, lambda = temp_mean)*rbinom(1,1,.5)
#     }
#   }
#
#   clone_fit <- emuFit_clone(Y = Y, X = X,
#                             method = "ML",
#                             reweight = FALSE)
#
#
#
#   non_cheat_fit <- emuFit(X = X, Y = Y,method = "ML",
#                           fast_fit = FALSE,
#                           fast_fit_cheat = FALSE,
#                           verbose = TRUE,
#                           reweight = FALSE,
#                           # weights= weights,
#                           tolerance = 1e-20)
#
#   expect_equal(non_cheat_fit$B, clone_fit$B,tolerance= 1e-3)
#
# })
#
# test_that("FL emuFit_clone returns same result as orig fit on sparse Y", {
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
#       Y[i,j] <- rpois(1, lambda = temp_mean)*rbinom(1,1,.5)
#     }
#   }
#
#   clone_fit <- emuFit_clone(Y = Y, X = X,
#                             method = "FL",
#                             reweight = FALSE)
#
#
#
#   non_cheat_fit <- emuFit(X = X, Y = Y,method = "FL",
#                           fast_fit = FALSE,
#                           fast_fit_cheat = FALSE,
#                           verbose = TRUE,
#                           reweight = FALSE,
#                           # weights= weights,
#                           tolerance = 1e-20)
#
#   expect_equal(non_cheat_fit$B, clone_fit$B,tolerance= 1e-4)
#
# })
#
#
# test_that("FL emuFit_clone returns same result as orig fit on sparse Y with
# z associated with B", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20))
#   z <- rnorm(40) + c(rep(3,20),rep(7,20))
#   b0 <- rnorm(10)
#   b1 <- 1:10
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 10, nrow = 40)
#
#   for(i in 1:40){
#     for(j in 1:10){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rpois(1, lambda = temp_mean)*rbinom(1,1,.5)
#     }
#   }
#
#   clone_fit <- emuFit_clone(Y = Y, X = X,
#                             method = "FL",
#                             reweight = FALSE)
#
#
#
#   non_cheat_fit <- emuFit(X = X, Y = Y,method = "FL",
#                           fast_fit = FALSE,
#                           fast_fit_cheat = FALSE,
#                           verbose = TRUE,
#                           reweight = FALSE,
#                           # weights= weights,
#                           tolerance = 1e-20)
#
#   expect_equal(non_cheat_fit$B, clone_fit$B,tolerance= 1e-2)
#
# })
#
# test_that("ML emuFit_clone returns same result as orig fit on sparse Y with
# z associated with B and more predictors", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20),rnorm(40),rnorm(40),rnorm(40))
#   z <- rnorm(40) + c(rep(3,20),rep(7,20))
#   b0 <- rnorm(10)
#   b1 <- 1:10
#   b <- rbind(b0,b1,b0,b1,b0)
#   Y <- matrix(NA,ncol = 10, nrow = 40)
#
#   for(i in 1:40){
#     for(j in 1:10){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rpois(1, lambda = temp_mean)*rbinom(1,1,.5)
#     }
#   }
#
#   clone_fit <- emuFit_clone(Y = Y, X = X,
#                             method = "ML",
#                             reweight = FALSE,
#                             maxit = 1000)
#
#
#
#   non_cheat_fit <- emuFit(X = X, Y = Y,method = "ML",
#                           fast_fit = FALSE,
#                           fast_fit_cheat = FALSE,
#                           verbose = TRUE,
#                           reweight = FALSE,
#                           # weights= weights,
#                           tolerance = 1e-20)
#
#   expect_equal(non_cheat_fit$B, clone_fit$B,tolerance= .2)
#
# })
#
# test_that("FL emuFit_clone returns same result as orig fit on sparse Y with
# z associated with B and there is separation", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20))
#   z <- rnorm(40) + c(rep(3,20),rep(7,20))
#   b0 <- rnorm(10)
#   b1 <- 1:10
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 10, nrow = 40)
#
#   for(i in 1:40){
#     for(j in 1:10){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rpois(1, lambda = temp_mean)*rbinom(1,1,.5 - j/21)
#     }
#   }
#
#   clone_fit <- emuFit_clone(Y = Y, X = X,
#                             method = "FL",
#                             reweight = FALSE)
#
#
#
#   non_cheat_fit <- emuFit(X = X, Y = Y,method = "FL",
#                           fast_fit = FALSE,
#                           fast_fit_cheat = FALSE,
#                           verbose = TRUE,
#                           reweight = FALSE,
#                           # weights= weights,
#                           tolerance = 1e-20)
#
#   expect_equal(non_cheat_fit$B, clone_fit$B,tolerance= 1e-2)
#
# })
#
#
#
