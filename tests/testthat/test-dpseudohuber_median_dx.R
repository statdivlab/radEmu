test_that("analytical deriv of pseudohuber_median wrt x matches numerical deriv", {

  set.seed(432)
  x <- rnorm(20)
  analytic_grad <- dpseudohuber_median_dx(x)
  numerical_grad <- numDeriv::grad(function(r) pseudohuber_median(x = r),
                                   x)

  expect_true(sum((analytic_grad - numerical_grad)^2)<1e-4)
})

test_that("analytical deriv of pseudohuber_median wrt x matches numerical deriv
for pseudohuber param d = 0.1", {

  set.seed(432)
  x <- rnorm(20)
  analytic_grad <- dpseudohuber_median_dx(x,d = 0.1)
  numerical_grad <- numDeriv::grad(function(r) pseudohuber_median(x = r,d = 0.1),
                                   x)

  expect_true(sum((analytic_grad - numerical_grad)^2)<1e-3)
})

test_that("analytical deriv of pseudohuber_median wrt x matches numerical deriv
for pseudohuber param d = 0.1", {

  set.seed(432)
  x <- c(1.944656e+01,
         -3.998268e+00,
         -3.025415e+00,
         -2.049766e+00,
         -1.031954e+00,
         -2.408555e-02,
         1.043626e+00 ,
         1.978452e+00 ,
         2.959190e+00,
         -2.484989e-09)
  analytic_grad <- dpseudohuber_median_dx(x,d = 0.1)
  numerical_grad <- numDeriv::grad(function(r) pseudohuber_median(x = r,d = 0.1),
                                   x)

  expect_true(sum((analytic_grad - numerical_grad)^2)<1e-3)
})
