test_that("Penalized estimation reduces to Haldane correction in saturated case", {
  Y <- matrix(1:4,nrow = 2)
  X <- cbind(1,c(0,1))
  penalized_fit <- emuFit_micro_penalized(X,
                                          Y,
                                          B = NULL,
                                          constraint_fn = function(x) x[2],
                                          maxit = 500,
                                          tolerance = 1e-5,
                                          verbose = FALSE)
  
  ml_fit <- emuFit_micro(X,Y + 0.5,
                         constraint_fn = function(x) x[2],
                         tolerance = 1e-5,
                         verbose = FALSE)
  
  expect_true(max(abs(penalized_fit$B- ml_fit)) <1e-4)
  expect_equal(penalized_fit$Y_augmented - Y, matrix(0.5, nrow = 2, ncol = 2))
  
  
  
  
  
})

test_that("Penalized estimation reduces to Haldane correction in saturated case with zeroes", {
  Y <- matrix(0:3,nrow = 2)
  X <- cbind(1,c(0,1))
  penalized_fit <- emuFit_micro_penalized(X,
                                          Y,
                                          B = NULL,
                                          constraint_fn = function(x) x[2],
                                          maxit = 500,
                                          tolerance = 1e-5,
                                          verbose = FALSE)
  
  ml_fit <- emuFit_micro(X,Y + 0.5,
                         constraint_fn = function(x) x[2],
                         tolerance = 1e-5,
                         verbose = FALSE)
  
  expect_true(max(abs(penalized_fit$B - ml_fit))<1e-4)
  expect_equal(penalized_fit$Y_augmented - Y, matrix(0.5, nrow = 2, ncol = 2))
})



test_that("PL fit to simple example returns reasonable values", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  z <- rnorm(40) +10
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
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }
  pl_fit <- emuFit_micro_penalized(X,
                                   Y,
                                   B = matrix(rnorm(20),nrow = 2),
                                   constraint_fn = function(x) mean(x),
                                   maxit = 200,
                                   tolerance = 1e-5,
                                   verbose= FALSE)
  
  # plot(b1-mean(b1),ml_fit[2,])
  # abline(a = 0,b = 1,lty =2,col = "red")
  
  expect_true(max(abs(pl_fit$B[2,] - (b1 - mean(b1))))<1e-2)
})


test_that("PL fit to simple example returns numerically identical results
regardless of whether we use computationally efficient augmentation or older
less efficient implementation (and that both substantially differ from MLE", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 5))
  z <- rnorm(10) +5
  J <- 10
  p <- 2
  n <- 10
  b0 <- rnorm(J)
  b1 <- seq(1,5,length.out = J)
  b1 <- b1 - mean(b1)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)
  
  for(i in 1:n){
    for(j in 1:J){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rnbinom(1, mu= temp_mean,size = 2)*rbinom(1,1,0.4)
    }
  }
  
  ml_fit <-  emuFit_micro(X,
                          Y,
                          B = NULL,
                          constraint_fn = function(x) mean(x),
                          maxit = 1000,
                          tolerance = 1e-3,
                          verbose= FALSE)
  pl_fit_new <- emuFit_micro_penalized(X,
                                       Y,
                                       B = NULL,
                                       constraint_fn = function(x) mean(x),
                                       maxit = 1000,
                                       tolerance = 1e-3,
                                       verbose= FALSE)
  #using old implementation should trigger a message explaining there's no
  #reason to do this except for testing
  suppressMessages(expect_message(
    pl_fit_old <- emuFit_micro_penalized(X,
                                         Y,
                                         B = NULL,
                                         constraint_fn = function(x) mean(x),
                                         maxit = 1000,
                                         tolerance = 1e-3,
                                         verbose= TRUE,
                                         use_legacy_augmentation = TRUE)))
  
  expect_true(max(abs(ml_fit - pl_fit_new$B))>1)
  expect_true(max(abs(ml_fit - pl_fit_old$B))>1)
  expect_equal(pl_fit_new,pl_fit_old)
})



test_that("PL fit to simple continuous covariate example returns numerically identical results
regardless of whether we use computationally efficient augmentation or older
less efficient implementation (and that both substantially differ from MLE", {
  set.seed(4323)
  X <- cbind(1,rnorm(10))
  z <- rnorm(10) +5
  J <- 10
  p <- 2
  n <- 10
  b0 <- rnorm(J)
  b1 <- seq(1,5,length.out = J)
  b1 <- b1 - mean(b1)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)
  
  for(i in 1:n){
    for(j in 1:J){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rnbinom(1, mu= temp_mean,size = 2)*rbinom(1,1,0.2)
    }
  }
  
  ml_fit <-  emuFit_micro(X,
                          Y,
                          B = NULL,
                          constraint_fn = function(x) mean(x),
                          maxit = 10000,
                          tolerance = 0.01,
                          verbose= FALSE)
  pl_fit_new <- emuFit_micro_penalized(X,
                                       Y,
                                       B = NULL,
                                       constraint_fn = function(x) mean(x),
                                       maxit = 10000,
                                       tolerance = 0.01,
                                       verbose= FALSE)
  #using old implementation should trigger a message explaining there's no
  #reason to do this except for testing
  suppressMessages(expect_message(
    pl_fit_old <- emuFit_micro_penalized(X,
                                         Y,
                                         B = NULL,
                                         constraint_fn = function(x) mean(x),
                                         maxit = 10000,
                                         tolerance = 0.01,
                                         verbose= TRUE,
                                         use_legacy_augmentation = TRUE)))
  
  expect_true(max(abs(ml_fit - pl_fit_new$B))>1)
  expect_true(max(abs(ml_fit - pl_fit_old$B))>1)
  expect_equal(pl_fit_new,pl_fit_old)
})




