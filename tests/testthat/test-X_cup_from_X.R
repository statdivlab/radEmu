test_that("X_cup_from_X_fast works and is faster than X_cup_from_X", {
  X <- cbind(1, rnorm(100), runif(100))
  J <- 200
  start1 <- proc.time()
  Xcup1 <- X_cup_from_X(X, J)
  end1 <- proc.time() - start1
  start2 <- proc.time()
  Xcup2 <- X_cup_from_X_fast(X, J)
  end2 <- proc.time() - start2
  
  expect_true(all.equal(Xcup1, Xcup2))
  expect_true(end2[3] < end1[3])
})
