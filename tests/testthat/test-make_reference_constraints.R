test_that("make_reference_constraints works", {
  # two columns in design matrix, reference taxon is taxon 5  
  list1 <- make_reference_constraints(p = 2, j = 5)

  # four columns in design matrix, reference taxon for all covariates is taxon 2
  list2 <- make_reference_constraints(p = 4, j = 2)

  # four columns in design matrix, reference taxon for covariates 1 and 2 is taxon 3 and
  # reference taxon for covariate 3 is taxon 4
  list3 <- make_reference_constraints(p = 4, j = c(3, 3, 4))
  
  expect_true(all.equal(length(list1), length(list2), length(list3), 2))
  expect_true(length(list1$constraint_list) == 2)
  expect_true(length(list2$constraint_list) == 4)
  expect_true(all.equal(list1$constraint_list[[1]](1:5), 
                        list2$constraint_list[[1]](1:5), 
                        list3$constraint_list[[1]](1:5)))
  expect_true(list1$constraint_list[[2]](2:6) == 6)
  expect_true(list3$constraint_list[[4]](11:15) == 14)
  expect_true(all.equal(list2$constraint_grad_list[[3]](1:4), c(0, 1, 0, 0)))
})
