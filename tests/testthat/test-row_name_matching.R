dat <- data.frame(cov1 = rep(c("A", "B", "C"), each = 6),
                  cov2 = rep(c("D", "E"), each = 9),
                  cov3 = rnorm(18),
                  cov4 = rnorm(18),
                  cov5 <- rep(c("G", "H", "I"), 6))

form <- ~ cov1 + cov2 + cov3 + cov4 + cov5
X.based <- model.matrix(form, dat)

Y <- matrix(rpois(18*6, 3), nrow = 18, ncol = 6)
colnames(Y) <- paste0("category_", 1:ncol(Y))
rownames(Y) <- paste0("sample_", 1:nrow(Y))

test_that("emuFit handles missing row names", {
  X1 <- matrix(X.based, nrow = 18)            # No row names
  
  expect_message(emuFit(Y = Y, X = X1, run_score_tests = FALSE), 
                 "Row names are missing from the covariate matrix X. We will assume the rows are in the same order as in the response matrix Y. You are responsible for ensuring the order of your observations is the same in both matrices.")
})

test_that("emuFit throws error on duplicate row names", {
  X2 <- X.based
  rownames(X2) <- rownames(Y)
  rownames(X2)[5] <- "sample_4" #Repeating one of the sample labels
  
  expect_error(emuFit(Y = Y, X = X2, run_score_tests = FALSE), 
               "Covariate matrix X has duplicated row names. Please ensure all row names are unique before refitting the model.")
})

test_that("emuFit subsets to common row names with warning", {
  X3 <- X.based
  rownames(X3) <- rownames(Y)
  X3 <- X3[c(1:4,7:14,16:18), , drop = FALSE]
  Y3 <- Y[c(1:2,5:18), , drop = FALSE]
  
  expect_warning(emuFit(Y = Y3, X = X3, run_score_tests = FALSE), 
                 regexp = "According to the rownames, there are observations that are missing either in the covariate matrix \\(X\\) and/or the response matrix \\(Y\\)\\. We will subset to common rows only, resulting in [0-9]+ samples\\.")
})

test_that("emuFit reorders rows of X when", {
  X4 <- X.based
  rownames(X4) <- rownames(Y)
  X4.p <- X4[c(1,(nrow(Y):2)), , drop = FALSE]

  expect_message(model.p <- emuFit(Y = Y, X = X4.p, run_score_tests = FALSE), 
                 "There is a different row ordering between the covariate matrix \\(X\\) and the response matrix \\(Y\\)\\. Covariate data will be reordered to match response data\\.")
  model.o <- emuFit(Y = Y, X = X4, run_score_tests = FALSE)
  expect_equal(model.o$coef$estimate, model.p$coef$estimate)
})

test_that("emuFit does not reorder rows of X when match_row_names is FALSE", {
  X5 <- X.based
  X5.p <- X.based
  rownames(X5) <- rownames(Y)
  rownames(X5.p) <- rownames(Y)[c(1,nrow(Y):2)]
  
  model.o <- emuFit(Y = Y, X = X5, run_score_tests = FALSE)
  model.p <- emuFit(Y = Y, X = X5.p, match_row_names = FALSE, run_score_tests = FALSE)
  
  expect_silent(model.o <- emuFit(Y = Y, X = X5, run_score_tests = FALSE))
  expect_silent(model.p <- emuFit(Y = Y, X = X5.p, match_row_names = FALSE, run_score_tests = FALSE))
  
  expect_equal(model.o$coef$estimate, model.p$coef$estimate)
})

test_that("emuFit stops when match_row_names is FALSE, but nrow does not coincide",{
  
  X6 <- X.based
  rownames(X6) <- rownames(Y)
  X6 <- X6[c(1,(nrow(Y):2)),]
  X6 <- X6[(1:16), , drop = FALSE]
  
  expect_error(emuFit(Y = Y, X = X6, match_row_names = FALSE, run_score_tests = FALSE), 
               "The number of rows does not match between the covariate matrix \\(X\\) and the response matrix \\(Y\\), and subsetting/matching by row name has been disabled\\. Please resolve this issue before refitting the model\\.")
})
