dat <- data.frame(cov1 = rep(c("A", "B", "C"), each = 6),
                  cov2 = rep(c("D", "E"), each = 9))
form1 <- ~cov1
form2 <- ~cov2
X1 <- model.matrix(form1, dat)
X2 <- model.matrix(form2, dat)
Y <- matrix(rpois(18*6, 3), nrow = 18, ncol = 6)
colnames(Y) <- paste0("category_", 1:ncol(Y))
Y0 <- Y
Y0[, c(1, 4)] <- 0
Y0[15, 1] <- 2
Y0[10, 4] <- 3

test_that("zero_comparison_check function works", {

  zero_comparison_res <- zero_comparison_check(X = X1, Y = Y0)
  expect_true(zero_comparison_res$zero_comparison[1])
  expect_true(zero_comparison_res$zero_comparison[10])
  expect_false(zero_comparison_res$zero_comparison[2])
  
})

test_that("zero_comparison column is added to coef when it should be", {
  
  # column doesn't exist with a 2 level covariate 
  emuRes1 <- emuFit(Y = Y0, X = X2, run_score_tests = FALSE, compute_cis = FALSE)
  expect_false("zero_comparison" %in% names(emuRes1$coef))
  
  # column doesn't exist when no zero-comparison parameters 
  emuRes2 <- emuFit(Y = Y, X = X1, run_score_tests = FALSE, compute_cis = FALSE)
  expect_false("zero_comparison" %in% names(emuRes2$coef))
  
  # column does exist when there are zero-comparison parameters
  emuRes3 <- emuFit(Y = Y0, X = X1, run_score_tests = FALSE, compute_cis = FALSE)
  expect_true("zero_comparison" %in% names(emuRes3$coef))
})

test_that("remove_zero_comparison_pvals argument works in different ways", {
  
  expect_error(emuFit(Y = Y0, X = X1, remove_zero_comparison_pvals = 15))
  
  emuRes1 <- emuFit(Y = Y0, X = X1, remove_zero_comparison_pvals = TRUE,
                    test_kj = data.frame(k = 2, j = 1))
  expect_true(is.na(emuRes1$coef$pval[1]))
  emuRes2 <- emuFit(Y = Y0, X = X1, remove_zero_comparison_pvals = TRUE,
                    test_kj = data.frame(k = 2, j = 1), use_fullmodel_info = TRUE,
                    return_both_score_pvals = TRUE)
  expect_true(is.na(emuRes2$coef$score_pval_full_info[1]) &  
                is.na(emuRes2$coef$score_pval_null_info[1]))
  emuRes3 <- emuFit(Y = Y0, X = X1, remove_zero_comparison_pvals = 0.5,
                    test_kj = data.frame(k = 2, j = 1:2))
  expect_true(is.na(emuRes3$coef$pval[1]) || emuRes3$coef$pval[1] > 0.5)
  expect_false(is.na(emuRes3$coef$pval[2]))
  emuRes4 <- emuFit(Y = Y0, X = X1, remove_zero_comparison_pvals = FALSE,
                    test_kj = data.frame(k = 2, j = 1))
  expect_true(emuRes4$coef$zero_comparison[1] & !is.na(emuRes4$coef$pval[1]))
  
})
