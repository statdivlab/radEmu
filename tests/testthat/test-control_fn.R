test_that("control_fn works as expected", {
  
  contr1 <- control_fn()
  contr2 <- control_fn(max_step = 1)
  expect_true(all.equal(contr1, contr2))
  
  contr3 <- control_fn(trackB = TRUE)
  expect_true(contr3$trackB)
  
  contr4 <- control_fn(contr3)
  expect_true(all.equal(contr3, contr4))
  
})
