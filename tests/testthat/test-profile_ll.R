test_that("profile_ll returns a number", {
  
  set.seed(88442)
  Y <- matrix(rexp(20),ncol = 5)
  X <- cbind(1,rep(0:1,2))
  B <- matrix(rnorm(10),nrow = 2)
  
  expect_true(is.numeric(profile_ll(X = X, Y = Y, B = B)))
})
