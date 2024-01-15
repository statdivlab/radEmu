test_that("sum_of_logs(x) numerically identical to log(sum(exp(x))", {
  set.seed(9043)
  x <- rnorm(100)
  expect_equal(sum_of_logs(x),log(sum(exp(x))))
})

test_that("sum_of_logs(c(x,y)) numerically identical to sum_of_logs(x,y) for x, y of length 1", {
  set.seed(934313)
  x <- rexp(1)
  y <- rexp(1)
  expect_equal(sum_of_logs(c(x,y)),sum_of_logs(x,y))
})

test_that("sum_of_logs(x) = x if x has length 1", {
  set.seed(93294)
  x <- rnorm(1)
  expect_equal(sum_of_logs(x),x)
})


test_that("sum_of_logs(x,y) with different length x and y returns error", {
  set.seed(900044)
  x <- 1:5
  y <- 1
  expect_error(sum_of_logs(x,y))
})


test_that("sum_of_logs(x) returns Inf if an element of x is Inf", {
  set.seed(9043)
  x <- c(2,Inf)
  expect_true(is.infinite(sum_of_logs(x)))
})

test_that("sum_of_logs(x) returns -Inf if  x = -Inf", {
  set.seed(9048853)
  x <- c(-Inf)
  expect_true(is.infinite(sum_of_logs(x)))
  expect_true(sum_of_logs(x)<0)
})

test_that("sum_of_logs(x) returns x_1 if x = c(x_1,-Inf)", {
  set.seed(2329043)
  x <- c(rnorm(1),-Inf)
  expect_equal(sum_of_logs(x),x[1])
})



test_that("sum_of_logs(c(x,-Inf) equal to sum_of_logs(x)", {
  set.seed(149043)
  x <-rnorm(100)
  expect_equal(sum_of_logs(c(x,-Inf)),sum_of_logs(x))
})

