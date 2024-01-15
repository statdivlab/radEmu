test_that("pseudohuber_center second deriv is equal to numerical derivative of first derivative", {

  set.seed(43234)
  x <- rnorm(10)

  ng <- numDeriv::grad(function(y) dpseudohuber_center_dx(y,d = 0.1)[1],
                 x)

  ag <- sapply(1:10,function(k) hess_pseudohuber_center(x,0.1,1,k))

  expect_equal(ng,ag,tolerance = 1e-5)
})
