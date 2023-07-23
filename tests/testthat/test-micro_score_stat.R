
test_that("Derivative matches analytical derivative", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  z <- rnorm(40) +3
  b0 <- rnorm(10)
  b1 <- seq(-1.5,1.5,length.out = 10)
  b1[5:6] <-  0
  b <- rbind(b0,b1)

  pvals <- numeric(100)
  score_stats <- numeric(100)
  Y <- matrix(NA,ncol = 10, nrow = 40)

    for(i in 1:40){
      for(j in 1:10){
        temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
        Y[i,j] <- rpois(1, lambda = temp_mean)
      }
    }

  n <- nrow(Y)
  J <- ncol(Y)
  p <- ncol(X)

  null_k <- 2
  null_j <- 5

  constrained_fit <- emuFit_micro_constrained(X = X,
                                              Y = Y,
                                              constrained_kj = c(null_k,null_j),
                                              maxit = maxit,
                                              constraint_fn = mean,
                                              tolerance = 0.01)

  score_stat <- micro_score_stat(Y = Y,
                                 X = X,
                                 B_fitted = constrained_fit$B,
                                 k = null_k,
                                 j = null_j
  )


  })
