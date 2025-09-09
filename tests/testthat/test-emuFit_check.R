set.seed(11)
J <- 6
n <- 12
X <- cbind(1,rnorm(n))
b0 <- rnorm(J)
b1 <- seq(1,5,length.out = J) -
  mean(seq(1,5,length.out = J))
b <- rbind(b0, b1)
Y <- radEmu:::simulate_data(n = n,
                            J = J,
                            X = X,
                            b0 = b0,
                            b1 = b1,
                            distn = "ZINB",
                            zinb_size = 2,
                            zinb_zero_prop = 0.2,
                            mean_z = 5)

#To ensure the messages about lack of row names do not show in the tests
rownames(X) <- paste0("Sample_",1:12)
rownames(Y) <- paste0("Sample_",1:12)

covariates <- data.frame(group = X[,2])

test_that("when X has a missing value, appropriate error hit", {
  X[2, 2] <- NA
  expect_error(emuFit(Y = Y, X = X), 
               "At least one value in your design matrix `X` is NA. Please remove any observations with missing values (NAs) or impute these missing values according to your analysis plan.",
               fixed = TRUE)
})

test_that("when X has a missing value, appropriate error hit", {
  covariates$group[2] <- NA
  suppressWarnings(expect_warning(
    emuFit(Y = Y, formula = ~ group, data = covariates, run_score_tests = FALSE),
    regexp = "Your data includes at least one NA value, and the `model.matrix()` function is automatically dropping all observations with any missing values. If you want to avoid this default behavior, please manually remove observations from your data with missing values or impute these missing values according to your analysis plan.",
    fixed = TRUE
  ))
})
