test_that("constraint_linearization is equal to constraint at point of linearization", {

constraint_fn <- pseudohuber_center
constraint_grad_fn <- dpseudohuber_center_dx
set.seed(4323)
x0 <- rnorm(10)
constraint_value <- constraint_fn(x0)

linearized_value <- constraint_linearization(x = x0,
                         x0 = x0,
                         constraint_fn_at_x0 = constraint_value,
                         constraint_grad_at_x0 = constraint_grad_fn(x0))

expect_equal(constraint_value,linearized_value)
})


test_that("constraint_linearization is equal to constraint in general when constraint is mean", {

  constraint_fn <- mean
  constraint_grad_fn <- (function(x) rep(1/length(x),length(x)))
  set.seed(4323)
  x0 <- rnorm(10)
  x <- rnorm(10)
  constraint_value_x0 <- constraint_fn(x0)

  linearized_value <- constraint_linearization(x = x,
                                               x0 = x0,
                                               constraint_fn_at_x0 = constraint_value_x0,
                                               constraint_grad_at_x0 = constraint_grad_fn(x0))

  constraint_value_x <- constraint_fn(x)
  expect_equal(constraint_value_x,linearized_value)
})


test_that("constraint_linearization is not equal to constraint in general
when constraint is nonlinear (pseudohuber), but that it approaches
correct value at point around which linear expansion is created", {

  constraint_fn <- pseudohuber_center
  constraint_grad_fn <- dpseudohuber_center_dx
  set.seed(4323)
  x0 <- rnorm(10)
  x <- rnorm(10)
  constraint_value_x0 <- constraint_fn(x0)

  linearized_value <- constraint_linearization(x = x,
                                               x0 = x0,
                                               constraint_fn_at_x0 = constraint_value_x0,
                                               constraint_grad_at_x0 = constraint_grad_fn(x0))

  constraint_value_x <- constraint_fn(x)
  expect_true(constraint_value_x != linearized_value)

  close_x <- x0 + 0.1*(x -x0)
  close_linearized_value <- constraint_linearization(x = close_x,
                                               x0 = x0,
                                               constraint_fn_at_x0 = constraint_value_x0,
                                               constraint_grad_at_x0 = constraint_grad_fn(x0))


  constraint_value_close_x <- constraint_fn(close_x)
  expect_true(sum(abs(close_linearized_value -   constraint_fn(close_x))) <
                sum(abs(linearized_value -    constraint_fn(x) )))


})
