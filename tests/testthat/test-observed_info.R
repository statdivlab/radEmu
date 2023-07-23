test_that("Observed info matches numerical hessian", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  z <- rnorm(40) +8
  J <- 10
  p <- 2
  n <- 40
  b0 <- rnorm(J)
  b1 <- 1:J
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
                         B = NULL,
                         constraint_fn = function(x) mean(x),
                         maxit = 200,
                         tolerance = 1e-5)

  obs_info <- observed_info(X,Y,ml_fit)

  numerical_hess <-
    numDeriv::hessian(function(x) -1*profile_ll_B_cup(X,Y,x),B_cup_from_B(ml_fit))


  expect_true(max(abs(obs_info/numerical_hess - 1)) < 1e-2)

})
