test_that("If I ask for silence, I get silence", {
  
  
  b1 <- 1:10
  b1 <- b1 - mean(b1)
  b1[5] <- pseudohuber_center(b1[-5],0.1)
  
  Y <- simulate_data(n=10, J=10, b0=rnorm(10), distn="Poisson", b1=b1, mean_count_before_ZI=500)
  
  expect_silent(full_fit <- emuFit(X = cbind(1, rep(c(0,1),each = 10/2)),
                                   Y = Y,
                                   B = NULL,
                                   tolerance = 1e-3, 
                                   verbose = FALSE))
  
  
})
