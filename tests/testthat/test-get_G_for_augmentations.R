test_that("get_G_for_augmentations_fast works and is faster than get_G_for_augmentations", {
  X <- cbind(1, rnorm(200), runif(200))
  J <- 200
  Xcup2 <- X_cup_from_X_fast(X, J)
  
  start1 <- proc.time()
  g1 <- get_G_for_augmentations(X, J, nrow(X), Xcup2)
  end1 <- proc.time() - start1
  
  start2 <- proc.time()
  g2 <- get_G_for_augmentations_fast(X, J, nrow(X), Xcup2)
  end2 <- proc.time() - start2
  
  expect_true(all.equal(g1, g2))
  expect_true(end2[3] < end1[3])
})
