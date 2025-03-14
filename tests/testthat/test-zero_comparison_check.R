dat <- data.frame(cov1 = rep(c("A", "B", "C"), each = 6),
                  cov2 = rep(c("D", "E"), each = 9),
                  cov3 = rnorm(18),
                  cov4 = rnorm(18),
                  cov5 <- rep(c("G", "H", "I"), 6))
form1 <- ~ cov1
form2 <- ~ cov2
form3 <- ~ cov1 + cov3
form4 <- ~ cov1 + cov3 + cov4
form5 <- ~ cov1 + cov2 + cov3 + cov4 + cov5
X1 <- model.matrix(form1, dat)
X1_base <- matrix(as.vector(X1), nrow = nrow(X1))
colnames(X1_base) <- colnames(X1)
X2 <- model.matrix(form2, dat)
X2_base <- matrix(as.vector(X2), nrow = nrow(X2))
colnames(X2_base) <- colnames(X2)
X3 <- model.matrix(form3, dat)
X3_base <- matrix(as.vector(X3), nrow = nrow(X3))
colnames(X3_base) <- colnames(X3)
X4 <- model.matrix(form4, dat)
X4_base <- matrix(as.vector(X4), nrow = nrow(X4))
colnames(X4_base) <- colnames(X4)
X5 <- model.matrix(form5, dat)
Y <- matrix(rpois(18*6, 3), nrow = 18, ncol = 6)
colnames(Y) <- paste0("category_", 1:ncol(Y))
Y0 <- Y
Y0[, c(1, 4)] <- 0
Y0[15, 1] <- 2
Y0[10, 4] <- 3

### check functionality when X does not have "assign" attribute from `model.matrix`
test_that("zero_comparison_check function works", {

  zero_comparison_res <- zero_comparison_check(X = X1_base, Y = Y0)
  expect_true(zero_comparison_res$zero_comparison[1])
  expect_true(zero_comparison_res$zero_comparison[10])
  expect_false(zero_comparison_res$zero_comparison[2])
  
})

test_that("zero_comparison column is added to coef when it should be", {
  
  # column doesn't exist with a 2 level covariate 
  emuRes1 <- emuFit(Y = Y0, X = X2_base, run_score_tests = FALSE, compute_cis = FALSE,
                    match_row_names = FALSE)
  expect_false("zero_comparison" %in% names(emuRes1$coef))
  
  # column doesn't exist when no zero-comparison parameters 
  emuRes2 <- emuFit(Y = Y, X = X1_base, run_score_tests = FALSE, compute_cis = FALSE,
                    match_row_names = FALSE)
  expect_false("zero_comparison" %in% names(emuRes2$coef))
  
  # column does exist when there are zero-comparison parameters
  emuRes3 <- emuFit(Y = Y0, X = X1_base, run_score_tests = FALSE, compute_cis = FALSE,
                    match_row_names = FALSE)
  expect_true("zero_comparison" %in% names(emuRes3$coef))
  
  # column does exist in the presence of continuous covariate 
  emuRes4 <- emuFit(Y = Y0, X = X4_base, run_score_tests = FALSE, compute_cis = FALSE,
                    match_row_names = FALSE)
  expect_true("zero_comparison" %in% names(emuRes4$coef))
  
})

test_that("remove_zero_comparison_pvals argument works in different ways", {
  
  expect_error(emuFit(Y = Y0, X = X1_base, remove_zero_comparison_pvals = 15,
                      match_row_names = FALSE))
  
  emuRes1 <- emuFit(Y = Y0, X = X1_base, remove_zero_comparison_pvals = TRUE,
                    test_kj = data.frame(k = 2, j = 1),
                    match_row_names = FALSE)
  expect_true(is.na(emuRes1$coef$pval[1]))
  emuRes2 <- emuFit(Y = Y0, X = X1_base, remove_zero_comparison_pvals = TRUE,
                    test_kj = data.frame(k = 2, j = 1), use_fullmodel_info = TRUE,
                    return_both_score_pvals = TRUE,
                    match_row_names = FALSE)
  expect_true(is.na(emuRes2$coef$score_pval_full_info[1]) &  
                is.na(emuRes2$coef$score_pval_null_info[1]))
  emuRes3 <- emuFit(Y = Y0, X = X1_base, remove_zero_comparison_pvals = 0.5,
                    test_kj = data.frame(k = 2, j = 1:2),
                    match_row_names = FALSE)
  expect_true(is.na(emuRes3$coef$pval[1]) || emuRes3$coef$pval[1] > 0.5)
  expect_false(is.na(emuRes3$coef$pval[2]))
  emuRes4 <- emuFit(Y = Y0, X = X1_base, remove_zero_comparison_pvals = FALSE,
                    test_kj = data.frame(k = 2, j = 1),
                    match_row_names = FALSE)
  expect_true(emuRes4$coef$zero_comparison[1] & !is.na(emuRes4$coef$pval[1]))
  
})

### check functionality when X does have "assign" attribute from `model.matrix`
test_that("zero_comparison column is added to coef when it should be, model.matrix X", {
  
  # column doesn't exist with a 2 level covariate 
  emuRes1 <- emuFit(Y = Y0, X = X2, run_score_tests = FALSE, compute_cis = FALSE,
                    match_row_names = FALSE)
  expect_false("zero_comparison" %in% names(emuRes1$coef))
  
  # column doesn't exist when no zero-comparison parameters 
  emuRes2 <- emuFit(Y = Y, X = X1, run_score_tests = FALSE, compute_cis = FALSE,
                    match_row_names = FALSE)
  expect_false("zero_comparison" %in% names(emuRes2$coef))
  
  # column does exist when there are zero-comparison parameters
  emuRes3 <- emuFit(Y = Y0, X = X1, run_score_tests = FALSE, compute_cis = FALSE,
                    match_row_names = FALSE)
  expect_true("zero_comparison" %in% names(emuRes3$coef))
  
  # column does exist in the presence of continuous covariate 
  emuRes4 <- emuFit(Y = Y0, X = X4, run_score_tests = FALSE, compute_cis = FALSE,
                    match_row_names = FALSE)
  expect_true("zero_comparison" %in% names(emuRes4$coef))
  
  # column does exist when there are multiple category covariates 
  emuRes5 <- emuFit(Y = Y0, X = X5, run_score_tests = FALSE, compute_cis = FALSE,
                    match_row_names = FALSE)
  expect_true("zero_comparison" %in% names(emuRes5$coef))
  expect_true(emuRes5$coef$zero_comparison[1])
  expect_false(emuRes5$coef$zero_comparison[13])
  expect_true(emuRes5$coef$zero_comparison[31])
  expect_false(emuRes5$coef$zero_comparison[32])
})
