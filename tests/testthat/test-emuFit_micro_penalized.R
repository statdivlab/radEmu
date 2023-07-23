test_that("Penalized estimation reduces to Haldane correction in saturated case", {
  Y <- matrix(1:4,nrow = 2)
  X <- cbind(1,c(0,1))
  penalized_fit <- emuFit_micro_penalized(X,
             Y,
             B = NULL,
             constraint_fn = function(x) x[2],
             maxit = 500,
             tolerance = 1e-5,
             collect_iterations = FALSE)

  ml_fit <- emuFit_micro(X,Y + 0.5,
                         constraint_fn = function(x) x[2],
                         tolerance = 1e-5)

  expect_true(max(abs(penalized_fit$B- ml_fit)) <1e-5)
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
                                          collect_iterations = FALSE)

  ml_fit <- emuFit_micro(X,Y + 0.5,
                         constraint_fn = function(x) x[2],
                         tolerance = 1e-5)

  expect_true(max(abs(penalized_fit$B - ml_fit))<1e-4)
  expect_equal(penalized_fit$Y_augmented - Y, matrix(0.5, nrow = 2, ncol = 2))
})



test_that("PL fit to simple example give reasonable output", {
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
                         tolerance = 1e-5)

  # plot(b1-mean(b1),ml_fit[2,])
  # abline(a = 0,b = 1,lty =2,col = "red")

  expect_true(max(abs(pl_fit$B[2,] - (b1 - mean(b1))))<1e-2)
})




test_that("PL fit to simple example give reasonable output under huber constraint", {
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
                                   constraint_fn = function(x) huber_center(x),
                                   maxit = 200,
                                   tolerance = 1e-5)

  # plot(b1-mean(b1),ml_fit[2,])
  # abline(a = 0,b = 1,lty =2,col = "red")

  expect_true(max(abs(pl_fit$B[2,] - (b1 - huber_center(b1))))< 1e-2)
})



test_that("PL fit to simple example give reasonable output under separation", {
  set.seed(43334)
  X <- cbind(1,rep(c(0,1),each = 20))
  z <- rep(-2,40)
  J <- 10
  p <- 2
  n <- 40
  b0 <- rnorm(J)
  b1 <- 2*seq(1,10,length.out = J)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = 40)

  for(i in 1:40){
    for(j in 1:J){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }
  rowSums(Y)
  colSums(Y[1:20,])
  pl_fit <- emuFit_micro_penalized(X,
                                   Y,
                                   B = matrix(rnorm(20),nrow = 2),
                                   constraint_fn = function(x) huber_center(x),
                                   maxit = 200,
                                   tolerance = 1e-2)


  # plot(b1 - huber_center(b1),pl_fit$B[2,])

 expect_true(mean((pl_fit$B[2,] - (b1 - huber_center(b1)))^2)<1)
})


test_that("PL fit to high-n low-count data suggests consistency", {
  set.seed(43334)
  n <- 200
  X <- cbind(1,rep(c(0,1),each = n/2))
  z <- rep(-1,n)
  J <- 10
  p <- 2
  b0 <- rnorm(J)
  b1 <- 2*seq(1,10,length.out = J)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)

  for(i in 1:n){
    for(j in 1:J){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }

  X <- X[rowSums(Y)>0,]
  Y <- Y[rowSums(Y)>0,]

  pl_fit <- emuFit_micro_penalized(X,
                                   Y,
                                   B = matrix(rnorm(2*J),nrow = 2),
                                   constraint_fn = function(x) huber_center(x),
                                   maxit = 200,
                                   tolerance = 1e-2)

#
#   plot(b1 - huber_center(b1),pl_fit$B[2,])
#   abline(a = 0,b = 1, lty = 2)

  expect_true(mean((pl_fit$B[2,] - (b1 - huber_center(b1)))^2)<.1)
})




