test_that("pseudohuber_center finds same center as numerical optimizer", {

  set.seed(4323)
  x <- rnorm(10)
  our_center <- pseudohuber_median(x)
  their_center <- optim(median(x),function(y) pseudohuber_loss(x- y,0.1),
                        method = "BFGS")$par

  expect_equal(our_center,their_center,tolerance = 1e-6)
})
