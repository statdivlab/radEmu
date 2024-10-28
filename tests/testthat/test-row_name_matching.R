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
  
  expect_message(emuFit(Y = Y, X = X1), 
                 "Row names are missing from the covariate matrix X. Assuming a one-to-one correspondence with the rows of the response matrix Y. Please double-check your data to confirm this correspondence.")
})

test_that("emuFit throws error on duplicate row names", {
  X2 <- X.based
  rownames(X2) <- rownames(Y)
  rownames(X2)[5] <- "sample_4" #Repeating one of the sample labels
  
  expect_error(emuFit(Y = Y, X = X2), 
               "Covariate matrix X has duplicated row names. Please ensure all row names are unique.")
})

test_that("emuFit subsets to common row names with warning", {
  X3 <- X.based
  rownames(X3) <- rownames(Y)
  X3 <- X3[c(1:4,7:14,16:18),]
  Y3 <- Y[c(1:2,5:18),]
  
  expect_warning(emuFit(Y = Y3, X = X3), 
                 regexp = "Row names differ between the covariate matrix \\(X\\) and the response matrix \\(Y\\)\\. Subsetting to common rows only, resulting in [0-9]+ samples\\.")
})

#------

test_that("emuFit reorders rows of X when", {
  X <- matrix(1:9, nrow = 3, dimnames = list(c("B", "C", "A"), NULL))
  Y <- matrix(9:1, nrow = 3, dimnames = list(c("A", "B", "C"), NULL))
  
  expect_warning(result <- MyFunc(X, Y, just_do_it = FALSE), 
                 "Row names do not match in order")
  expect_equal(rownames(result$X), rownames(result$Y))  # Ensure rows are reordered
})

test_that("MyFunc does not reorder rows of X when just_do_it is TRUE", {
  X <- matrix(1:9, nrow = 3, dimnames = list(c("B", "C", "A"), NULL))
  Y <- matrix(9:1, nrow = 3, dimnames = list(c("A", "B", "C"), NULL))
  
  result <- MyFunc(X, Y, just_do_it = TRUE)
  expect_equal(rownames(result$X), c("B", "C", "A"))   # Original order maintained
  expect_equal(rownames(result$Y), c("A", "B", "C"))
})
