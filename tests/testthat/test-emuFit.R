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

