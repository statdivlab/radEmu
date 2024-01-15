test_that("Conversion between B and B_cup recovers same matrix", {
  set.seed(0)
  B <- matrix(rnorm(40),ncol = 10)
  B_cup = B_cup_from_B(B)
  B_prime = B_from_B_cup(B_cup, J = ncol(B), p = nrow(B))
  B_cup_prime = B_cup_from_B(B_prime)

  expect_equal(B, B_prime)
  expect_equal(B_cup, B_cup_prime)
})
