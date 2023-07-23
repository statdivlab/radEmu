test_that("emuFit works on simple example",{
    set.seed(4323)
    X <- cbind(1,rep(c(0,1),each = 20))
    z <- rnorm(40) +5
    J <- 10
    p <- 2
    n <- 40
    b0 <- rnorm(J)
    b1 <- -.1*seq(1,10,length.out = J)
    b <- rbind(b0,b1)
    Y <- matrix(NA,ncol = J, nrow = 40)

    for(i in 1:40){
      for(j in 1:J){
        temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
        Y[i,j] <- rpois(1, lambda = temp_mean)
      }
    }

    data <- data.frame(covariate = X[,2])
    fitted_emu <- emuFit(Y,~covariate,data)

    expect_true(is.data.frame(fitted_emu))

})


test_that("emuFit works on simple example with wider variation in counts by taxon",{
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  z <- rnorm(40) +5
  J <- 10
  p <- 2
  n <- 40
  b0 <- rnorm(J)
  b1 <- seq(1,10,length.out = J)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = 40)

  for(i in 1:40){
    for(j in 1:J){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rnbinom(1, mu = temp_mean,size = 1.2)#rpois(1, lambda = temp_mean)
    }
  }

  data <- data.frame(covariate = X[,2])
  fitted_emu <- emuFit(Y,~covariate,data)
  expect_true(is.data.frame(fitted_emu))
})


test_that("emuFit works on example with more taxa",{
  set.seed(4323)
  n <- 6
  X <- cbind(1,rep(c(0,1),each = n/2))
  z <- rnorm(n) +5
  J <- 30
  p <- 2

  b0 <- rnorm(J)
  b1 <- seq(1,10,length.out = J)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)

  for(i in 1:n){
    for(j in 1:J){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rnbinom(1, mu = temp_mean,size = 1.2)#rpois(1, lambda = temp_mean)
    }
  }

  data <- data.frame(covariate = X[,2])
  fitted_emu <- emuFit(Y,~covariate,data)
  expect_true(is.data.frame(fitted_emu))
})


# test_that("Formulas work", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20))
#   covariate_data <- data.frame("group" = rep(c(0,1),each = 20))
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
#   ml_fit <- emuFit(formula_rhs = ~group,
#                    Y = Y,
#                    covariate_data = covariate_data,
#                    method = "ML")
#
#
#
#   expect_true(max(abs(ml_fit$B[2,] - seq(-4.5,4.5,1)))<.02)
# })
#
#
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
#   ml_fit <- emuFit(X = X, Y = Y,method = "ML",
#                                      tolerance = 0.01)
#
#
#
#   expect_true(max(abs(ml_fit$B[2,] - seq(-4.5,4.5,1)))<.02)
# })
#
# test_that("ML fit to simple example give reasonable output under null hypothesis", {
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
#   ml_fit <- emuFit(X = X, Y = Y,method = "ML",
#                    tolerance = 0.01,
#                    constraint_fn = pseudohuber_center,
#                    optim_only = FALSE)
#   B <- ml_fit$B
#   B[2,] <- B[2,] - pseudohuber_center(B[2,-5])
#   B[2,5] <- 0
#   ml_restricted_fit <-
#     emuFit(X = X,
#            Y = Y,
#            method = "ML",
#            tolerance = 0.01,
#            maxit = 1e5,
#            B = B,
#            constraint_fn = pseudohuber_center,
#            constrained_k = 2, #currently supports only 1 constrained k
#            constrained_j = 5,
#            verbose = TRUE,
#            optim_only = FALSE)
#
#   ml_restricted_fit$B
#
#
#   expect_true(max(abs(ml_fit$B[2,] - seq(-4.5,4.5,1)))<.02)
# })
#
# test_that("Fast fit ML using a `cheat` works", {
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
#   cheat_fit <- emuFit(X = X, Y = Y,method = "ML",
#                    fast_fit = TRUE,
#                    fast_fit_cheat = TRUE,
#                    verbose = TRUE,
#                    tolerance = 0.01)
#
#
#   expect_true(max(abs(cheat_fit$B[2,] - seq(-4.5,4.5,1)))<.02)
# })
#
# test_that("Fast fit ML with weights using a `cheat` works", {
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
#   weights <- matrix(rgamma(40,shape= 1/10,rate = 1),ncol = 1)%*%matrix(1,ncol = 10,
#                                                 nrow = 1)
#
#   cheat_fit <- emuFit(X = X, Y = Y,method = "ML",
#                       fast_fit = TRUE,
#                       fast_fit_cheat = TRUE,
#                       verbose = TRUE,
#                       reweight = FALSE,
#                       weights = weights,
#                       tolerance = 0.01)
#
#   non_cheat_fit <- emuFit(X = X, Y = Y,method = "ML",
#                           fast_fit = FALSE,
#                           fast_fit_cheat = FALSE,
#                           verbose = TRUE,
#                           reweight = FALSE,
#                           weights = weights,
#                           tolerance = 0.01)
#
#   expect_equal(cheat_fit$B,non_cheat_fit$B,tolerance = 1e-2)
# })
#
#
#
# test_that("ML with B_[,J] = 0 gives same result
# up to shift as median-constrained fit (both with NB data)", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20))
#   z <- rnorm(40) +7
#   b0 <- seq(-10,10,length.out = 10)
#   b1 <- 1:10
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 10, nrow = 40)
#
#   for(i in 1:40){
#     for(j in 1:10){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rnbinom(1, mu = temp_mean,size = 0.5)
#     }
#   }
#   corner_fit <- emuFit(X = X, Y = Y,method = "ML",
#                      fast_fit = FALSE,
#                      fast_fit_cheat = FALSE,
#                      constraint_fn = (function(x) x[10]),
#                      verbose = TRUE,
#                      tolerance = 0.01)
#
#   median_fit <- emuFit(X = X, Y = Y,method = "ML",
#                          fast_fit = FALSE,
#                          fast_fit_cheat = FALSE,
#                          verbose = TRUE,
#                          tolerance = 0.01)
#
#   corner_fit$B[2,] <- corner_fit$B[2,]-median(corner_fit$B[2,])
#   corner_fit$B[1,] <- corner_fit$B[1,]-median(corner_fit$B[1,])
#   expect_true(max(abs(corner_fit$B[2,] - seq(-4.5,4.5,1)))<2)
#   expect_equal(corner_fit$B,median_fit$B,tolerance = 1e-4)
#   expect_equal(corner_fit$weights,median_fit$weights,tolerance = 1e-4)
# })
#
# test_that("ML with B_[,J] = 0 gives same result
# up to shift as median-constrained fit (both with NB data and weights)", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20))
#   z <- rnorm(40) +7
#   b0 <- seq(-10,10,length.out = 10)
#   b1 <- 1:10
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 10, nrow = 40)
#
#   for(i in 1:40){
#     for(j in 1:10){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rnbinom(1, mu = temp_mean,size = 5)
#     }
#   }
#   weights <- matrix(rgamma(40,shape = 1/100, scale = 1),ncol = 1)%*%
#     matrix(1, nrow = 1, ncol = 10)
#   corner_fit <- emuFit(X = X, Y = Y,method = "ML",
#                        fast_fit = FALSE,
#                        fast_fit_cheat = FALSE,
#                        constraint_fn = (function(x) x[10]),
#                        verbose = TRUE,
#                        optim_only= FALSE,
#                        weights = weights,
#                        tolerance = 0.01)
#
#   median_fit <- emuFit(X = X, Y = Y,method = "ML",
#                        fast_fit = FALSE,
#                        fast_fit_cheat = FALSE,
#                        verbose = TRUE,
#                        optim_only = FALSE,
#                        weights = weights,
#                        tolerance = 0.01)
#
#
#
#
#   corner_fit$B[2,] <- corner_fit$B[2,]-median(corner_fit$B[2,])
#   corner_fit$B[1,] <- corner_fit$B[1,]-median(corner_fit$B[1,])
#   expect_true(max(abs(corner_fit$B[2,] - seq(-4.5,4.5,1)))<20)
#   expect_equal(corner_fit$B[2,],median_fit$B[2,],tolerance = 1e-4)
#   expect_equal(corner_fit$weights,median_fit$weights,tolerance = 1e-4)
# })
#
# test_that("emuFit_clone and emuFit return same results with NB data", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20))
#   z <- rnorm(40) +7
#   b0 <- seq(-5,5,length.out = 10)
#   b1 <- 1:10 - 5.5
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 10, nrow = 40)
#
#   for(i in 1:40){
#     for(j in 1:10){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rnbinom(1, mu = temp_mean,size = 5)
#     }
#   }
#
#   median_fit <- emuFit(X = X, Y = Y,method = "ML",
#                        fast_fit = FALSE,
#                        fast_fit_cheat = FALSE,
#                        verbose = TRUE,
#                        optim_only = FALSE,
#                        tolerance = 0.01)
#
#   median_clone_fit <- emuFit_clone(Y = Y,
#                                    X = X)
#
#   expect_equal(median_clone_fit$B[2,],median_fit$B[2,],tolerance = 0.25)
# })
#
# test_that("Fast fit FL *not* using a `cheat` works with NB data", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20))
#   z <- rnorm(40) +7
#   b0 <- seq(-10,10,length.out = 10)
#   b1 <- 1:10
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 10, nrow = 40)
#
#   for(i in 1:40){
#     for(j in 1:10){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rnbinom(1, mu = temp_mean,size = 0.5)
#     }
#   }
#   fast_fit <- emuFit(X = X, Y = Y,method = "FL",
#                      fast_fit = TRUE,
#                      fast_fit_cheat = TRUE,
#                      fast_fit_anchor = 10,
#                      verbose = TRUE,
#                      tolerance = 0.01)
#
#   non_fast_fit <- emuFit(X = X, Y = Y,method = "FL",
#                          fast_fit = FALSE,
#                          fast_fit_cheat = FALSE,
#                          verbose = TRUE,
#                          tolerance = 0.01)
#
#
#   expect_true(max(abs(fast_fit$B[2,] - seq(-4.5,4.5,1)))<2)
#   expect_equal(fast_fit$B,non_fast_fit$B,tolerance = 1e-4)
# })
#
#
# test_that("Fast fit reweighted ML using a `cheat` works with NB data", {
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
#   cheat_fit <- emuFit(X = X, Y = Y,method = "ML",
#                       fast_fit = TRUE,
#                       fast_fit_cheat = TRUE,
#                       verbose = TRUE,
#                       reweight = TRUE,
#                       tolerance = 0.01)
#
#   non_cheat_fit <- emuFit(X = X, Y = Y,method = "ML",
#                           fast_fit = FALSE,
#                           fast_fit_cheat = FALSE,
#                           verbose = TRUE,
#                           reweight = TRUE,
#                           tolerance = 0.01)
#
#
#   expect_true(max(abs(cheat_fit$B[2,] - seq(-4.5,4.5,1)))<3)
#   expect_equal(cheat_fit$B,non_cheat_fit$B,tolerance = 1e-4)
#   expect_equal(cheat_fit$weights,non_cheat_fit$weights,tolerance = 1e-4)
# })
#
#
# test_that("Fast fit weighted FL using a `cheat` works", {
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
#   weights <- matrix(rgamma(40,shape = 1/100,scale = 1), ncol = 1) %*%
#     matrix(1,ncol = 10,nrow = 1)
#
#   cheat_fit <- emuFit(X = X, Y = Y,method = "FL",
#                       fast_fit = TRUE,
#                       fast_fit_cheat = TRUE,
#                       verbose = TRUE,
#                       reweight = FALSE,
#                       weights = weights,
#                       tolerance = 1e-2)
#
#   non_cheat_fit <- emuFit(X = X, Y = Y,method = "FL",
#                           fast_fit = FALSE,
#                           fast_fit_cheat = FALSE,
#                           verbose = TRUE,
#                           reweight = FALSE,
#                           weights = weights,
#                           tolerance = 0.01)
#
#   # plot(non_cheat_fit$B[2,]-cheat_fit$B[2,])
#
#   expect_true(max(abs(cheat_fit$B[2,] - seq(-4.5,4.5,1)))<.1)
#   expect_equal(cheat_fit$B,non_cheat_fit$B,tolerance = 1e-2)
#   expect_equal(cheat_fit$weights,non_cheat_fit$weights,tolerance = 1e-3)
# })
#
# test_that("Fast fit FL using a `cheat` works", {
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
#   cheat_fit <- emuFit(X = X, Y = Y,method = "FL",
#                       fast_fit = TRUE,
#                       fast_fit_cheat = TRUE,
#                       verbose = TRUE,
#                       tolerance = 0.01)
#
#   non_cheat_fit <- emuFit(X = X, Y = Y,method = "FL",
#                           fast_fit = FALSE,
#                           fast_fit_cheat = FALSE,
#                           verbose = TRUE,
#                           tolerance = 0.01)
#
#
#   expect_true(max(abs(cheat_fit$B[2,] - seq(-4.5,4.5,1)))<.02)
#   expect_equal(cheat_fit$B,non_cheat_fit$B,tolerance = 1e-4)
# })
#
# test_that("Fast fit ML using a `cheat` works when some entries of Y are 0", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20))
#   z <- rnorm(40)
#   z[1] <- 0
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
#   cheat_fit <- emuFit(X = X, Y = Y,method = "ML",
#                       fast_fit = TRUE,
#                       fast_fit_cheat = TRUE,
#                       verbose = TRUE,
#                       tolerance = 0.01)
#
#   non_cheat_fit <-
#     emuFit(X = X, Y = Y,method = "ML",
#            fast_fit = FALSE,
#            fast_fit_cheat = FALSE,
#            verbose = TRUE,
#            tolerance = 0.01)
#
#
#   expect_true(max(abs(cheat_fit$B[2,] - seq(-4.5,4.5,1)))<.5)
#   expect_equal(cheat_fit$B,non_cheat_fit$B,tolerance = 1e-4)
# })
#
#
#
# test_that("ML fit to with multiple covariates give reasonable output", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20),rnorm(20))
#   z <- rnorm(40) +8
#   b0 <- rnorm(10)
#   b1 <- 1:10
#   b2 <- 1:10
#   b <- rbind(b0,b1,b2)
#   Y <- matrix(NA,ncol = 10, nrow = 40)
#
#   for(i in 1:40){
#     for(j in 1:10){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rpois(1, lambda = temp_mean)
#     }
#   }
#   ml_fit <- emuFit(X = X, Y = Y,method = "ML",
#                    tolerance = 0.01)
#
#   expect_true(max(abs(ml_fit$B[2,] - seq(-4.5,4.5,1)))<.05)
#   expect_true(max(abs(ml_fit$B[3,] - seq(-4.5,4.5,1)))<.01)
# })
#
# test_that("reweighted ML fit to simple example give reasonable output", {
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
#                                      reweight = TRUE,
#                                      tolerance = 0.01)
#   expect_true(max(abs(ml_fit$B[2,] - seq(-4.5,4.5,1)))<.02)
# })
#
#
# test_that("FL fit to simple example give reasonable output", {
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
#   fl_fit <- fitted_model_fl <- emuFit(X = X, Y = Y,method = "FL")
#   expect_true(max(abs(fl_fit$B[2,] - seq(-4.5,4.5,1)))<.02)
# })
#
# test_that("Reweighted FL fit to simple example give reasonable output", {
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
#   fl_fit <- fitted_model_fl <- emuFit(X = X, Y = Y,method = "FL",reweight= TRUE)
#   expect_true(max(abs(fl_fit$B[2,] - seq(-4.5,4.5,1)))<.02)
# })
#
# #
# # test_that("FL fit to simple example give reasonable output in
# # larger example", {
# #   set.seed(4323)
# #   X <- cbind(1,rep(c(0,1),each = 200))
# #   z <- rnorm(400) +8
# #   b0 <- rnorm(300)
# #   b1 <- seq(1,10,length.out = 300)
# #   b <- rbind(b0,b1)
# #   Y <- matrix(NA,ncol = 300, nrow = 400)
# #
# #   for(i in 1:400){
# #     for(j in 1:300){
# #       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
# #       Y[i,j] <- rpois(1, lambda = temp_mean)
# #     }
# #   }
# #   fl_fit <- fitted_model_fl <- emuFit(X = X, Y = Y,method = "FL")
# #   expect_true(max(abs(fl_fit$B[2,] - (b1 - median(b1))))<0.01)
# # })
#
# test_that("FL fit to very simple example give reasonable output", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 5))
#   z <- rnorm(10) +8
#   b0 <- rnorm(2)
#   b1 <- 1:2
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 2, nrow = 10)
#
#   for(i in 1:10){
#     for(j in 1:2){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rpois(1, lambda = temp_mean)
#     }
#   }
#   fl_fit <- fitted_model_fl <- emuFit(X = X, Y = Y,method = "FL")
#   expect_true(max(abs(fl_fit$B[2,] -c(-.5,.5)))<0.01)
# })
#
# test_that("FL fit same with different identifiability constraints
# (up to location shift)", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 5))
#   z <- rnorm(10) +8
#   b0 <- rnorm(2)
#   b1 <- 1:2
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 2, nrow = 10)
#
#   for(i in 1:10){
#     for(j in 1:2){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rpois(1, lambda = temp_mean)
#     }
#   }
#   fl_fit_med <- emuFit(X = X, Y = Y,method = "FL")
#   fl_fit_j <- emuFit(X = X, Y = Y,method = "FL",
#                      constraint_fn = (function(x)x[1]))
#   expect_equal(fl_fit_j$B[,1] - fl_fit_j$B[,2] ,
#        fl_fit_med$B[,1] - fl_fit_med$B[,2],
#        tolerance = 0.01)
# })
#
#
# test_that("FL fit to simple example same up to location shift when
# different identifiability constraints used", {
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
#   fl_fit_med <- emuFit(X = X, Y = Y,method = "FL")
#   fl_fit_j <- emuFit(X = X, Y = Y,method = "FL",
#                      constraint_fn = function(x){ x[1]})
#   expect_equal(fl_fit_med$B[2,] -fl_fit_med$B[2,1],
#                fl_fit_j$B[2,] -fl_fit_j$B[2,1],
#                tolerance = 0.01)
# })
#
# test_that("Reweighted ML fit to simple example give reasonable output", {
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
#                                      tolerance = 0.01,
#                                      reweight = TRUE)
#   expect_true(max(abs(ml_fit$B[2,] - seq(-4.5,4.5,1)))<.02)
# })
#
#
# test_that("Reweighted FL fit to simple example give reasonable output", {
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
#   fl_fit <- fitted_model_fl <- emuFit(X = X, Y = Y,method = "FL",
#                                       reweight = TRUE)
#   expect_true(max(abs(fl_fit$B[2,] - seq(-4.5,4.5,1)))<.02)
# })
#
#
# test_that("FL fit in case with separation yields finite estimates", {
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
#   Y[1:20,10] <- 0
#   fl_fit <- fitted_model_fl <- emuFit(X = X, Y = Y,method = "FL")
#   expect_true(max(abs(fl_fit$B[2,1:9] - seq(-4.5,3.5,1)))<.03)
#   # expect_true(fl_fit$B[2,10]<17)
# })
#
# test_that("Directly computing hat matrix yields same data augmentation
# as Cholesky decomposition of information.", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20))
#   p <- ncol(X)
#   n <- nrow(X)
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
#   Y[,10] <- 0
#   J <- ncol(Y)
#
#   alt_augmentations <- emuFit_clone(X = X,
#                                     Y = Y,
#                                     method = "FL")
#
#   augmentations <- fitted_model_fl <- emuFit(X = X, Y = Y,method = "FL",
#                                              B = alt_augmentations$B,
#                                       test_firth = TRUE)
#
#
#   alt_augmentations <- alt_augmentations$Y_augmented - alt_augmentations$Y
#
#  expect_equal(augmentations$augmentations_naive,
#               augmentations$augmentations_chol,tolerance = 1e-5)
#
#  expect_equal(augmentations$augmentations_naive,
#              as.numeric(Y_to_Y_tilde(alt_augmentations)),
#               tolerance = 1e-5)
#
#  expect_equal(sum(alt_augmentations),(n+(J - 1)*p)/2)
#
#
# })
#
# test_that("Directly computing hat matrix yields same data augmentation
# as Cholesky decomposition of information, but, like, actually compute
# the hat matrix, David.", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20))
#   z <- rnorm(40) +8
#   b0 <- rnorm(10)
#   b1 <- 1:10
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 10, nrow = 40)
#   J <- ncol(Y)
#   n <- nrow(Y)
#   p <- ncol(X)
#
#   for(i in 1:40){
#     for(j in 1:10){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rpois(1, lambda = temp_mean)
#     }
#   }
#   Y[,10] <- 0
#   alt_augmentations <- emuFit_clone(X = X,
#                                     Y = Y,
#                                     method = "FL")
#
#   weights <- as.numeric(alt_augmentations$weights)
#
#   augmentations <- fitted_model_fl <- emuFit(X = X, Y = Y,method = "FL",
#                                              B = alt_augmentations$B,
#                                              test_firth = TRUE)
#
#   B <- alt_augmentations$B
#   z <- alt_augmentations$z
#   alt_augmentations <- alt_augmentations$Y_augmented - alt_augmentations$Y
#
#   expect_equal(augmentations$augmentations_naive,
#                augmentations$augmentations_chol,tolerance = 1e-5)
#
#   X_tilde <- X_to_X_tilde(X,J)
#   Y_tilde <- Y_to_Y_tilde(Y)
#   S <- Matrix::sparseMatrix(i = 1:(n*J),
#                             j = rep(1:n,each = J),
#                             x = rep(1, n*J))
#   D_tilde <- cbind(X_tilde,S)
#
#   rownames(D_tilde) <- 1:(n*J)
#   B_tilde <- B_to_B_tilde(B)
#   theta <- rbind(B_tilde,Matrix::Matrix(z,ncol = 1))
#   X_tilde_repar <- X_tilde
#   X_tilde_repar <- X_tilde_repar[,-(p*J - (0:(p-1)))]
#   D_tilde_repar <- cbind(X_tilde_repar,S)
#
#   B_tilde <- B_to_B_tilde(B)
#   theta <- rbind(B_tilde,Matrix::Matrix(z,ncol = 1))
#   W <- Matrix::Diagonal(x = weights*as.numeric(exp((D_tilde%*%theta))))
#
#   theta <- rbind(B_tilde,Matrix::Matrix(z,ncol = 1))
#   W_half <- Matrix::Diagonal(x =
#                                sqrt(as.numeric(weights*exp((D_tilde%*%theta)))))
#   # info <- Matrix::crossprod(D_tilde_repar,W_half)%*%W_half%*%D_tilde_repar
#   augmentations_direct <- (W_half%*%D_tilde_repar)%*%qr.solve(
#     as.matrix(Matrix::crossprod(W_half%*%D_tilde_repar)),
#       as.matrix(Matrix::crossprod(D_tilde_repar,W_half)),tol = 1e-20)
#
#   augmentations_direct <- 0.5*diag(as.matrix(augmentations_direct))
#
#   expect_equal(augmentations_direct, augmentations$augmentations_chol,
#                tolerance = .5)
#
# })
#
# test_that("FL fit to simple sparse data yields finite solutions
# and otherwise is similar to ML fit", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 5))
#   z <- rnorm(10) +3
#   b0 <- rnorm(10)
#   b1 <- 1:10
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 10, nrow = 10)
#
#   for(i in 1:10){
#     for(j in 1:10){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rpois(1, lambda = temp_mean)
#     }
#   }
#   Y[6:10,4:6] <- 0
#   ml_fit <- emuFit(X = X, Y = Y,method = "ML",
#                                      tolerance = 0.01)
#   fl_fit <- emuFit(X = X, Y = Y,method = "FL",
#                    # constraint = function(x) x[10],
#                    tolerance = 0.01)
#   expect_equal(ml_fit$B[,-c(4:6)],fl_fit$B[,-c(4:6)],tolerance = 0.1)
#   expect_true(sum(abs(ml_fit$B[,4:6]))/sum(abs(fl_fit$B[,4:6]))>2)
# })
#
# test_that("FL fit in case with separation actually maximizes PL", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20))
#   z <- rnorm(40) +8
#   b0 <- rnorm(10)
#   b1 <- 1:10
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 10, nrow = 40)
#   n <- 40
#   J <- 10
#   p <- 2
#
#   for(i in 1:40){
#     for(j in 1:10){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rpois(1, lambda = temp_mean)
#     }
#   }
#   Y[1:20,10] <- 0
#   fl_fit <- fitted_model_fl <- emuFit(X = X, Y = Y,method = "FL",
#                                       tolerance = 1e-3)
#
#   X_tilde <- X_to_X_tilde(X,J)
#   S <-  S <- Matrix::sparseMatrix(i = 1:(n*J),
#                                   j = rep(1:n,each = J),
#                                   x = rep(1, n*J))
#   D_tilde <- cbind(X_tilde,S)
#   pl_fn <- function(x){
#     B_temp <- fl_fit$B
#     B_temp[2,10] <-     B_temp[2,10] + x
#   B_tilde <- B_to_B_tilde(B_temp)
#   theta <- rbind(B_tilde,Matrix::Matrix(z,ncol = 1))
#
#   return(log_likelihood_wide(Y,
#                              fl_fit$weights,
#                              X,
#                              B_temp,
#                              fl_fit$z) +
#     calculate_firth_penalty(D_tilde = D_tilde,
#                             W = Matrix::Diagonal(x =
#                                                    as.numeric(
#                                                      exp(D_tilde%*%theta))),
#                             n_skip = 2))
#   }
#
#   ds <- seq(-.5,.5,by = .01)
#   pls <- sapply(ds,pl_fn)
#
#   expect_equal(ds[which.max(pls)],0)
#
#
# })
#
