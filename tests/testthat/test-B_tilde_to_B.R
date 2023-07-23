# test_that("B_to_B_tilde and B_tilde_to_B are inverses",{
#   set.seed(0)
#   p <- 5
#   J <- 17
#   B <- matrix(rnorm(p*J),
#               nrow = p,
#               ncol = J)
#   B_tilde <- B_to_B_tilde(B)
#   B_prime <- B_tilde_to_B(B_tilde,J = J, p = p)
#   B_tilde_prime <- B_to_B_tilde(B_prime)
#   expect_equal(B,B_prime)
#   expect_equal(B_tilde,B_tilde_prime)
# })
