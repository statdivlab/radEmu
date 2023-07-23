test_that("analytical hessian closely approximates expected information", {


  Y <- cbind(c(4,6),c(2,3))
  Y_long <- rbind(matrix(Y[1,],ncol = 1),
                  matrix(Y[2,],ncol = 1))
  X <- cbind(1,c(0,1))
  X_tilde <- X_cup_from_X(X,2)
  n <- 2
  J <- 2
  p <- 2

  X_tilde_J <- X_tilde
  for(i in 1:n){
    X_tilde_J[(i - 1)*J + J,] <- 0*X_tilde_J[(i - 1)*J + J,]
  }
  X_tilde_J <- X_tilde_J[,1:(p*(J - 1))]

  B <- emuFit_micro(X,Y,constraint_fn = function(x) x[2],tolerance = 1e-5)

  tau <- rowSums(Y)

  numerical_hessian <-
    numDeriv::hessian(function(x) log_likelihood_repar(Y_long,
                                                       beta_tilde_J = matrix(x[1:2],ncol = 1),
                                                       tau = x[3:4],
                                                       X_tilde_J,
                                                       J,
                                                       n,
                                                       p),
                      x = c(B[,1],tau))


  analytical_info <- expected_info_repar(beta_tilde_J = B[,1,drop = FALSE],
                                  tau = tau,
                                  X_tilde_J = X_tilde_J,
                                  J = 2,
                                  n = 2,
                                  p = 2,
                                  collect_byproducts = TRUE)$expected_info

  expect_true(max(abs(analytical_info + numerical_hessian)) <1e-5)


})
