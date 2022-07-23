test_that("Conversion of Y to and from long format yields same matrix", {
  set.seed(0)
  J <- 10
  Y <- matrix(rnorm(200),ncol = J)
  Y_tilde <- Y_to_Y_tilde(Y)
  Y_prime <- Y_tilde_to_Y(Y_tilde,J)
  expect_equal(Y,Y_prime)
})

test_that("Conversion of Y_tilde to and from wide format yields same matrix", {
  set.seed(0)
  J <- 10
  Y_tilde<- matrix(rnorm(200),ncol = 1)
  Y <- Y_tilde_to_Y(Y_tilde,J)
  Y_tilde_prime <- Y_to_Y_tilde(Y)
  expect_equal(Y_tilde,Y_tilde_prime)
})
