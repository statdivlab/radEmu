
test_that("ML fit to simple example give reasonable output", {
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
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }
  ml_fit <- emuFit_micro(X,
             Y,
             constraint_fn = function(x) mean(x),
             maxit = 200,
             tolerance = 1e-2)

  # plot(b1-mean(b1),ml_fit[2,])
  # abline(a = 0,b = 1,lty =2,col = "red")

  expect_true(max(abs(ml_fit[2,] - (b1 - mean(b1))))<.1)
})


test_that("ML fit to multiple regressors, large n, and excess-Poisson variance gives reasonable output", {
  set.seed(4323)
  n <- 800
  X <- cbind(1,rep(c(0,1),each = n/2),rnorm(n),rnorm(n),rnorm(n))
  z <- rnorm(n) +5
  J <- 10
  p <- 2

  b0 <- rnorm(J)
  b1 <- seq(1,10,length.out = J)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)

  for(i in 1:n){
    for(j in 1:J){
      temp_mean <- exp(X[i,1:2,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rnbinom(1,mu= temp_mean,size = 0.25)#rpois(1, lambda = temp_mean)
    }
  }
  ml_fit <- emuFit_micro(X,
                         Y,
                         constraint_fn = function(x) mean(x),
                         maxit = 200,
                         tolerance = 1e-2)

  # plot(b1-mean(b1),ml_fit[2,])
  # abline(a = 0,b = 1,lty =2,col = "red")

  expect_true(max(abs(ml_fit[2,] - (b1 - mean(b1))))<.5)
})



test_that("ML fit to multiple regressors, larger n, truly high J, and excess-Poisson variance gives reasonable output", {
  set.seed(4323)
  n <- 800
  X <- cbind(1,rep(c(0,1),each = n/2),rnorm(n),rnorm(n),rnorm(n))
  z <- rnorm(n) +5
  J <- 1000
  p <- 2

  b0 <- rnorm(J)
  b1 <- seq(1,10,length.out = J)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)

  for(i in 1:n){
    for(j in 1:J){
      temp_mean <- exp(X[i,1:2,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rnbinom(1,mu= temp_mean,size = 0.25)#rpois(1, lambda = temp_mean)
    }
  }
  init_fit <- emuFit_micro_penalized(X,
                         Y,
                         maxit = 500,
                         tolerance = 1e2,
                         max_step = 0.1,
                         verbose = TRUE)

  # plot(b1-mean(b1),ml_fit[2,])
  # abline(a = 0,b = 1,lty =2,col = "red")

  expect_true(max(abs(ml_fit[2,] - (b1 - mean(b1))))<.5)
})


test_that("ML fit to multiple regressors, moderate n, large J, and excess-Poisson variance gives reasonable output", {
  set.seed(4323)
  n <- 40
  X <- cbind(1,rep(c(0,1),each = n/2),rnorm(n),rnorm(n),rnorm(n))
  z <- rnorm(n)
  J <- 100
  p <- 2

  b0 <- rnorm(J)
  b1 <- seq(1,10,length.out = J)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)

  for(i in 1:n){
    for(j in 1:J){
      temp_mean <- exp(X[i,1:2,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rnbinom(1,mu= temp_mean,size = 0.25)#rpois(1, lambda = temp_mean)
    }
  }
  ml_fit <- emuFit_micro(X,
                         Y,
                         constraint_fn = function(x) mean(x),
                         maxit = 200,
                         tolerance = 1e-1,
                         max_step = 0.5)

  # plot(b1-mean(b1),ml_fit[2,])
  # abline(a = 0,b = 1, lty = 2, col = "red")

  b1_hat <- ml_fit[2,]
  expect_true(abs(lm(b1_hat~b1)$coef[2] - 1) < 0.2)
  # abline(a = 0,b = 1,lty =2,col = "red")

})

# t(ml_fit) %>% as.data.frame() %>% mutate(j = 1:J) %>%
#   pivot_longer(-j) %>%
#   ggplot() +
#   geom_point(aes(x = j,y = value, color = name)) +
#   geom_line(aes(x = j, y  = value, group = name, color = name)) +
#   theme_bw()



test_that("ML fit to simple example give reasonable output with J > n", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  z <- rnorm(40) +8
  b0 <- rnorm(80)
  b1 <- seq(-5,5,length.out = 80)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = 80, nrow = 40)

  for(i in 1:40){
    for(j in 1:80){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }
  ml_fit <- emuFit_micro(X,
                         Y,
                         B = NULL,
                         constraint_fn = function(x) mean(x),
                         maxit = 500,
                         tolerance = 0.01)





  expect_true(max(abs(ml_fit[2,] - b1))<.1)

  # X_cup <- X_cup_from_X(X,10)
  # B_cup <- B_cup_from_B(ml_fit)



})



test_that("ML fit to simple example give reasonable output with J >> n", {
  set.seed(4323)
  n <- 10
  X <- cbind(1,rep(c(0,1),each = n/2))
  z <- rnorm(n) +8
  J <- 1000
  b0 <- rnorm(J)
  b1 <- seq(-5,5,length.out = J)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)

  for(i in 1:n){
    for(j in 1:J){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }


  ml_fit <- emuFit_micro(X,
                         Y,
                         constraint_fn = function(x) mean(x),
                         maxit = 500,
                         tolerance = 0.01)





  expect_true(max(abs(ml_fit[2,] - b1))<.5)


  # plot(b1,ml_fit[2,],pch = ".")


})





