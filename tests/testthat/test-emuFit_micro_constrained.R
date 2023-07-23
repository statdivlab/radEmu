test_that("Constrained fit returns correct values when null is true", {
  set.seed(4323)
  J = 40
  p = 2
  n <- 20
  X <- cbind(1,rep(c(0,1),each = n/2))
  z <- rnorm(n) +8
  # b0 <- c(-1,0,1)
  # b1 <- c(-1,0,1)
  b0 <- rnorm(J)
  b1 <- seq(-3,3,length.out = J)
  b1[2] <- pseudohuber_center(b1[-2])
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)


  for(i in 1:n){
    for(j in 1:J){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
      # Y[i,j] <- rnbinom(1,size = 1.25,mu = temp_mean)
    }
  }

  unconstrained_fit <-
    emuFit_micro(X =X, Y = Y)


  constrained_fit <-
    emuFit_micro_constrained(X,
           Y,
           B = unconstrained_fit,
           constraint_fn = mean,# pseudohuber_center,
           constraint_grad_fn = function(x) rep(1/length(x),length(x)),#pseudohuber_center_dx,
           k_constr = 2,
           j_constr = 2,
           maxit = 500,
           info_reg = 0.5,
           tolerance = 0.01)


  expect_equal(constrained_fit[2,],b[2,] - mean(b[2,]),tolerance= 0.01)

})

test_that("Constrained fit returns correct values when null is true", {
  set.seed(4323)
  J = 40
  p = 2
  n <- 20
  X <- cbind(1,rep(c(0,1),each = n/2))
  z <- rnorm(n) +8
  # b0 <- c(-1,0,1)
  # b1 <- c(-1,0,1)
  b0 <- rnorm(J)
  b1 <- seq(-3,3,length.out = J)
  b1[2] <- pseudohuber_center(b1[-2])
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)


  for(i in 1:n){
    for(j in 1:J){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
      # Y[i,j] <- rnbinom(1,size = 5,mu = temp_mean)
    }
  }

  unconstrained_fit <-
    emuFit_micro(X =X, Y = Y)


  constrained_fit <-
    emuFit_micro_constrained(X,
                             Y,
                             B = unconstrained_fit,
                             constraint_fn = pseudohuber_center,
                             constraint_grad_fn = dpseudohuber_center_dx,
                             k_constr = 2,
                             j_constr = 2,
                             maxit = 500,
                             info_reg = 0.5,
                             tolerance = 0.01)


  expect_equal(constrained_fit[2,],b[2,] - pseudohuber_center(b[2,]),tolerance= 0.01)

})


test_that("Constrained fit returns correct values when null is true and irrelevant covariates are present", {
  set.seed(4323)
  J = 40
  p = 2
  n <- 20
  X <- cbind(1,rep(c(0,1),each = n/2))
  z <- rnorm(n) +8
  # b0 <- c(-1,0,1)
  # b1 <- c(-1,0,1)
  b0 <- rnorm(J)
  b1 <- seq(-3,3,length.out = J)
  b1[2] <- pseudohuber_center(b1[-2])
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)


  for(i in 1:n){
    for(j in 1:J){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
      # Y[i,j] <- rnbinom(1,size = 5,mu = temp_mean)
    }
  }

  X <- cbind(X,rnorm(n),rnorm(n))

  unconstrained_fit <-
    emuFit_micro(X =X, Y = Y)


  constrained_fit <-
    emuFit_micro_constrained(X,
                             Y,
                             B = unconstrained_fit,
                             constraint_fn = pseudohuber_center,
                             constraint_grad_fn = dpseudohuber_center_dx,
                             k_constr = 2,
                             j_constr = 2,
                             maxit = 500,
                             info_reg = 0.5,
                             tolerance = 0.01)


  expect_equal(constrained_fit[2,],b[2,] - pseudohuber_center(b[2,]),tolerance= 0.01)

})

