test_that("analytical and numerical derivatives of profile ll align at mle and at
another value of B", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  z <- rnorm(40) +8
  b0 <- rnorm(10)
  b1 <- 1:10
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = 10, nrow = 40)
  J <- 10
  p <-  2

  for(i in 1:40){
    for(j in 1:10){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }
  ml_fit <- emuFit_micro(X,
                         Y,
                         B = NULL,
                         constraint_fn = function(x) mean(x),
                         maxit = 500,
                         tolerance = 1e-3)

  B <- ml_fit
  B_cup <- B_cup_from_B(B)
  analytical <- dpll_dB_cup(X,Y,B)
  for_numderiv <- function(B_cup){
    profile_ll(X,
               Y,
               B_from_B_cup(matrix(B_cup,ncol = 1),J,p))}

  numerical <- numDeriv::grad(for_numderiv,as.numeric(B_cup))

  expect_true(max(abs(analytical-numerical))<0.1)

  B <- B + rnorm(20)

  B_cup <- B_cup_from_B(B)

  analytical <- dpll_dB_cup(X,Y,B)

  numerical <- numDeriv::grad(for_numderiv,as.numeric(B_cup))


  expect_true(max(abs(analytical/numerical) -1)<1e-5)


})
