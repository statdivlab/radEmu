test_that("pseudohuber_center finds same center as numerical optimizer", {

  set.seed(4323)
  x <- rnorm(10)
  our_center <- pseudohuber_center(x)
  their_center <- optim(median(x),function(y) pseudohuber(x- y,1),
                        method = "BFGS")$par

  expect_equal(our_center,their_center,tolerance = 1e-6)
})
