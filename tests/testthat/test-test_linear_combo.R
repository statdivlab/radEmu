set.seed(11)
J <- 7
n <- 12
X <- cbind(1,rnorm(n), rnorm(n))
b0 <- rnorm(J)
b1 <- seq(1,5,length.out = J) -
  mean(seq(1,5,length.out = J))
b2 <- seq(1,7,length.out = J) -
  mean(seq(1,7,length.out = J))
b <- rbind(b0, b1, b2)
Y <- radEmu:::simulate_data(n = n,
                            J = J,
                            X = X,
                            B = b, 
                            distn = "ZINB",
                            zinb_size = 2,
                            zinb_zero_prop = 0.2,
                            mean_z = 5)

rownames(X) <- paste0("Sample_",1:12)
rownames(Y) <- paste0("Sample_",1:12)
covariates <- data.frame(cov1 = X[, 2], cov2 = X[, 3])
emuRes <- emuFit(X = X, Y = Y, run_score_tests = FALSE, tolerance = 0.01, return_wald_p = TRUE)

test_that("test-linear-combo gives same result as Wald test when linear combo is c(0, 1)", {
  lin_combo_test <- test_linear_combo(fitted_model = emuRes, linear_combo = c(0, 1, 0), j = 1)
  expect_true(all.equal(unname(emuRes$coef[1, c("estimate", "lower", "upper", "wald_p")]),
                        unname(lin_combo_test[c("linear_combo_estimate", "ci_lower", 
                                                   "ci_upper", "wald_pval")])))
})

test_that("test-linear-combo can run a linear combo with multiple values", {
  lin_combo_test1 <- test_linear_combo(fitted_model = emuRes, linear_combo = c(0, 1, -1), j = 1)
  lin_combo_test4 <- test_linear_combo(fitted_model = emuRes, linear_combo = c(0, 1, -1), j = 4)
  lin_combo_test14 <- test_linear_combo(fitted_model = emuRes, linear_combo = c(0, 1, -1), j = c(1, 4))
  expect_true(lin_combo_test1$wald_pval < lin_combo_test4$wald_pval)
  expect_equal(lin_combo_test1$linear_combo_estimate, emuRes$B[2, 1] - emuRes$B[3, 1])
  all.equal(lin_combo_test14[1, ], lin_combo_test1)
})

# EDIT THIS ONE!!

# test reparameterization - scale
test_that("test-linear-combo with reparameterization", {
  X1 <- X
  X1[, 2] <- 100*X[, 2]
  emuRes1 <- emuFit(X = X1, Y = Y, run_score_tests = FALSE, tolerance = 0.01, return_wald_p = TRUE)
  lin_combo_test <- test_linear_combo(fitted_model = emuRes, linear_combo = c(0, 1/100, 0))
})

test_that("linear combo test controls t1e for Poisson data, n = 100", {
  
  skip(message = "Skipping test t1e of Wald linear combos because this takes a long time")
  
  nsim <- 100
  set.seed(11)
  J <- 11
  n <- 100
  X <- cbind(1,rnorm(n), rnorm(n))
  b0 <- rnorm(J)
  b1 <- seq(1,5,length.out = J) -
    mean(seq(1,5,length.out = J))
  b2 <- seq(1,7,length.out = J) -
    mean(seq(1,7,length.out = J))
  b2[4] <- -.8
  b2[8] <- .8
  b <- rbind(b0, b1, b2)
  
  res <- data.frame(sim = 1:nsim, j1 = NA, j4 = NA, j5 = NA, j6 = NA, j8 = NA)
  for (sim in 1:nsim) {
    print(sim)
    Y <- radEmu:::simulate_data(n = n, J = J,
                                X = X,
                                B = b, 
                                distn = "Poisson",
                                zinb_size = 2,
                                zinb_zero_prop = 0.2,
                                mean_z = 5)
    emuRes <- emuFit(X = X, Y = Y, run_score_tests = FALSE, tolerance = 0.001, return_wald_p = TRUE,
                     match_row_names = FALSE)
    lin_combo_test <- test_linear_combo(fitted_model = emuRes, linear_combo = c(0, 1, -1), j = c(1, 4, 5, 6, 8))
    res[sim, c(2:6)] <- lin_combo_test$wald_pval
  }
  
  expect_false(mean(res$j1 <= 0.05) <= 0.05) # b1 = -2, b2 = -3
  expect_true(mean(res$j4 <= 0.05) <= 0.05) # b1 = b2 = -0.8
  expect_false(mean(res$j5 <= 0.05) <= 0.05) # b1 = -0.4, b2 = -0.6
  expect_true(mean(res$j6 <= 0.05) <= 0.05) # b1 = b2 = 0
  expect_true(mean(res$j8 <= 0.05) <= 0.05) # b1 = b2 = 0.8
  
})

# EDIT THIS ONE!!
# test reparameterization - location
set.seed(11)
J <- 7
n <- 120
cov <- rep(c("A", "B", "C"), each = 40)
b0 <- rnorm(J)
b0 <- b0 - radEmu:::pseudohuber_center(b0, 0.1)
b1 <- rnorm(J)
b1 <- b1 - radEmu:::pseudohuber_center(b1, 0.1)
b2 <- rnorm(J)
b2 <- b2 - radEmu:::pseudohuber_center(b2, 0.1)
b <- rbind(b0, b1, b2)
X1 <- model.matrix(~ cov, data = data.frame(cov))
X2 <- model.matrix(~ cov, data = data.frame(cov))
Y <- radEmu:::simulate_data(n = n,
                            J = J,
                            X = X1,
                            B = b, 
                            distn = "ZINB",
                            zinb_size = 2,
                            zinb_zero_prop = 0.2,
                            mean_z = 5)
emuRes1 <- emuFit(X = X1, Y = Y, run_score_tests = FALSE, tolerance = 0.01, return_wald_p = TRUE,
                  match_row_names = FALSE)
emuRes2 <- emuFit(X = X2, Y = Y, run_score_tests = FALSE, tolerance = 0.01, return_wald_p = TRUE,
                  match_row_names = FALSE)


test_that("test-linear-combo with reparameterization", {
  lin_combo_test <- test_linear_combo(fitted_model = emuRes1, linear_combo = c(1, 1, 0))
})

