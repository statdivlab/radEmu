test_that("Test analytical information's closeness to numerical hessian at
MLE in two-group model", {
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
                         tolerance = 1e-4)

  B <- ml_fit
  B_cup <- B_cup_from_B(B)
  analytical_info <- f_info(Y,B_cup,X,X_cup)
  analytical <- dpll_dB_cup(X,Y,B)
  for_numderiv <- function(B_cup){
    profile_ll(X,
               Y,
               B_from_B_cup(matrix(B_cup,ncol = 1),J,p))}

  numerical <- -1*numDeriv::hessian(for_numderiv,as.numeric(B_cup))

 expect_true(max(abs(as.matrix(analytical_info)/numerical - 1))< 1e-2)
 expect_true(max(abs(asinh(as.matrix(analytical_info)) - asinh(numerical )))<1e-2)
})
