test_that("Pseudodeterminant approximates determinant in full rank case", {
  set.seed(432359)
  A <- matrix(rnorm(100),ncol = 10)

  direct_log_det <- Matrix::determinant(A)$modulus
  pseudo_det <- log_pseudodet(A,mat_rank = 10)

  expect_equal(direct_log_det,pseudo_det,tolerance = 1e-3)
})

test_that("Pseudodeterminant approximates determinant of full rank submatrix
in rank deficient case", {
  set.seed(5e3)
  A <- matrix(rnorm(100),ncol = 10)
  direct_log_det <- Matrix::determinant(A)$modulus

  A <- cbind(A,0)
  A <- rbind(A,0)


  pseudo_det <- log_pseudodet(A,mat_rank = 10)

  expect_equal(direct_log_det,pseudo_det,tolerance = 1e-3)
})
