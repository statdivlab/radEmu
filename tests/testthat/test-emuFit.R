test_that("ML fit to simple example give reasonable output", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  z <- rnorm(40) +8
  b0 <- rnorm(10)
  b1 <- 1:10
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = 10, nrow = 40)

  for(i in 1:40){
    for(j in 1:10){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }
  ml_fit <- emuFit(X = X, Y = Y,method = "ML",
                                     tolerance = 0.01)

  lin_pred <- ml_fit$X%*%ml_fit$B + matrix(ml_fit$z,nrow = 40,ncol = 1)%*%
    matrix(1,nrow = 1, ncol = 10)

  prefitted <- exp(lin_pred)


  expect_true(max(abs(ml_fit$B[2,] - seq(-4.5,4.5,1)))<.01)
})

test_that("reweighted ML fit to simple example give reasonable output", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  z <- rnorm(40) +8
  b0 <- rnorm(10)
  b1 <- 1:10
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = 10, nrow = 40)

  for(i in 1:40){
    for(j in 1:10){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }
  ml_fit <- fitted_model_ml<- emuFit(X = X, Y = Y,method = "ML",
                                     reweight = TRUE,
                                     tolerance = 0.01)
  expect_true(max(abs(ml_fit$B[2,] - seq(-4.5,4.5,1)))<.01)
})


test_that("FL fit to simple example give reasonable output", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  z <- rnorm(40) +8
  b0 <- rnorm(10)
  b1 <- 1:10
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = 10, nrow = 40)

  for(i in 1:40){
    for(j in 1:10){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }
  fl_fit <- fitted_model_fl <- emuFit(X = X, Y = Y,method = "FL")
  expect_true(max(abs(fl_fit$B[2,] - seq(-4.5,4.5,1)))<.01)
})

test_that("Reweighted FL fit to simple example give reasonable output", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  z <- rnorm(40) +8
  b0 <- rnorm(10)
  b1 <- 1:10
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = 10, nrow = 40)

  for(i in 1:40){
    for(j in 1:10){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }
  fl_fit <- fitted_model_fl <- emuFit(X = X, Y = Y,method = "FL",reweight= TRUE)
  expect_true(max(abs(fl_fit$B[2,] - seq(-4.5,4.5,1)))<.01)
})

#
# test_that("FL fit to simple example give reasonable output in
# larger example", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 200))
#   z <- rnorm(400) +8
#   b0 <- rnorm(300)
#   b1 <- seq(1,10,length.out = 300)
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 300, nrow = 400)
#
#   for(i in 1:400){
#     for(j in 1:300){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rpois(1, lambda = temp_mean)
#     }
#   }
#   fl_fit <- fitted_model_fl <- emuFit(X = X, Y = Y,method = "FL")
#   expect_true(max(abs(fl_fit$B[2,] - (b1 - median(b1))))<0.01)
# })

test_that("FL fit to very simple example give reasonable output", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 5))
  z <- rnorm(10) +8
  b0 <- rnorm(2)
  b1 <- 1:2
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = 2, nrow = 10)

  for(i in 1:10){
    for(j in 1:2){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }
  fl_fit <- fitted_model_fl <- emuFit(X = X, Y = Y,method = "FL")
  expect_true(max(abs(fl_fit$B[2,] -c(-.5,.5)))<0.01)
})

test_that("FL fit same with different identifiability constraints
(up to location shift)", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 5))
  z <- rnorm(10) +8
  b0 <- rnorm(2)
  b1 <- 1:2
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = 2, nrow = 10)

  for(i in 1:10){
    for(j in 1:2){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }
  fl_fit_med <- emuFit(X = X, Y = Y,method = "FL")
  fl_fit_j <- emuFit(X = X, Y = Y,method = "FL",
                     constraint_fn = (function(x)x[1]))
  expect_equal(fl_fit_j$B[,1] - fl_fit_j$B[,2] ,
       fl_fit_med$B[,1] - fl_fit_med$B[,2],
       tolerance = 0.01)
})


test_that("FL fit to simple example same up to location shift when
different identifiability constraints used", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  z <- rnorm(40) +8
  b0 <- rnorm(10)
  b1 <- 1:10
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = 10, nrow = 40)

  for(i in 1:40){
    for(j in 1:10){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }
  fl_fit_med <- emuFit(X = X, Y = Y,method = "FL")
  fl_fit_j <- emuFit(X = X, Y = Y,method = "FL",
                     constraint_fn = function(x){ x[1]})
  expect_equal(fl_fit_med$B[2,] -fl_fit_med$B[2,1],
               fl_fit_j$B[2,] -fl_fit_j$B[2,1],
               tolerance = 0.01)
})

test_that("Reweighted ML fit to simple example give reasonable output", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  z <- rnorm(40) +8
  b0 <- rnorm(10)
  b1 <- 1:10
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = 10, nrow = 40)

  for(i in 1:40){
    for(j in 1:10){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }
  ml_fit <- fitted_model_ml<- emuFit(X = X, Y = Y,method = "ML",
                                     tolerance = 0.01,
                                     reweight = TRUE)
  expect_true(max(abs(ml_fit$B[2,] - seq(-4.5,4.5,1)))<.01)
})


test_that("Reweighted FL fit to simple example give reasonable output", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  z <- rnorm(40) +8
  b0 <- rnorm(10)
  b1 <- 1:10
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = 10, nrow = 40)

  for(i in 1:40){
    for(j in 1:10){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }
  fl_fit <- fitted_model_fl <- emuFit(X = X, Y = Y,method = "FL",
                                      reweight = TRUE)
  expect_true(max(abs(fl_fit$B[2,] - seq(-4.5,4.5,1)))<.01)
})


test_that("FL fit in case with separation yields finite estimates", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  z <- rnorm(40) +8
  b0 <- rnorm(10)
  b1 <- 1:10
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = 10, nrow = 40)

  for(i in 1:40){
    for(j in 1:10){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }
  Y[1:20,10] <- 0
  fl_fit <- fitted_model_fl <- emuFit(X = X, Y = Y,method = "FL")
  expect_true(max(abs(fl_fit$B[2,1:9] - seq(-4.5,3.5,1)))<.01)
  expect_true(fl_fit$B[2,10]<15)
})

test_that("Directly computing hat matrix yields same data augmentation
as Cholesky decomposition of information.", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  z <- rnorm(40) +8
  b0 <- rnorm(10)
  b1 <- 1:10
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = 10, nrow = 40)

  for(i in 1:40){
    for(j in 1:10){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }
  Y[,10] <- 0
  augmentations <- fitted_model_fl <- emuFit(X = X, Y = Y,method = "FL",
                                      test_firth = TRUE)

 expect_equal(augmentations$augmentations_naive,
              augmentations$augmentations_chol)
})

test_that("FL fit to simple sparse data yields finite solutions
and otherwise is similar to ML fit", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 5))
  z <- rnorm(10) +3
  b0 <- rnorm(10)
  b1 <- 1:10
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = 10, nrow = 10)

  for(i in 1:10){
    for(j in 1:10){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }
  Y[6:10,4:6] <- 0
  ml_fit <- emuFit(X = X, Y = Y,method = "ML",
                                     tolerance = 0.01)
  fl_fit <- emuFit(X = X, Y = Y,method = "FL",
                   tolerance = 0.01)
  expect_equal(ml_fit$B[,-c(4:6)],fl_fit$B[,-c(4:6)],tolerance = 0.1)
  expect_true(sum(abs(ml_fit$B[,4:6]))/sum(abs(fl_fit$B[,4:6]))>2)
})
