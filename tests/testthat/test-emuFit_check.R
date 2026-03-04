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

covariates <- data.frame(group = X[,2],
                         cat = rep(c("A", "B", "C"), each = 4))
other_covariates <- model.matrix(~cat, covariates)

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

test_that("when k in test_kj is a string it works", {
  expect_silent(emuFit(Y = Y, formula = ~ group, data = covariates, 
                       test_kj = data.frame(j = 2, k = "group")))
  res1 <- emuFit(Y = Y, formula = ~ cat, data = covariates, 
                       test_kj = data.frame(j = 2, k = "cat"))
  expect_true(sum(!(is.na(res1$coef$score_stat))) == 2)
  res2 <- emuFit(Y = Y, formula = ~ cat, data = covariates, 
                test_kj = data.frame(j = 2, k = "catB"))
  expect_true(all.equal(res1$coef[2, ], res2$coef[2, ]))
  res3 <- emuFit(Y = Y, formula = ~ group + cat, data = covariates, 
                 test_kj = data.frame(j = 2, k = c(2, "cat")))
  expect_true(sum(!(is.na(res3$coef$score_stat))) == 3)
})

  