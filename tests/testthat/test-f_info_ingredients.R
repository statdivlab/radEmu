
test_that("Fisher info constructed from X_cup and W matches info constructed by summing over i directly.", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  z <- rnorm(40) +8
  J <- 5
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
                         B = NULL,
                         constraint_fn = function(x) mean(x),
                         maxit = 200,
                         tolerance = 1e-1)

  B_cup <- B_cup_from_B(ml_fit)
  X_cup <- X_cup_from_X(X,J)
  direct_sum_info <- f_info(Y,B_cup, B = ml_fit,X,X_cup)
  W <- f_info_ingredients(Y,B_cup,B = ml_fit,X,X_cup)

  W_info <- Matrix::t(X_cup)%*%W%*%X_cup

  expect_true(max(abs(W_info -direct_sum_info)) < 1e-5)



})

