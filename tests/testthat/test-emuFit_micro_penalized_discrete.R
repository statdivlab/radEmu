test_that("closed-form discrete augmentations equal get_augmentations at arbitrary B", {
  set.seed(1)
  n <- 15
  J <- 6
  g <- sample(rep(1:3, c(3, 5, 7)))
  X <- cbind(1, g == 2, g == 3)
  Y <- matrix(rpois(n * J, 5), n, J)
  B <- matrix(rnorm(3 * J), 3, J)
  G <- get_G_for_augmentations_fast(X, J, n, X_cup_from_X_fast(X, J))
  expect_lt(max(abs(get_augmentations(X = X, G = G, Y = Y, B = B) -
                      get_augmentations_discrete(X = X, Y = Y, B = B))), 1e-10)
})

test_that("closed-form discrete penalized fit matches iterative fit, binary design with separation", {
  set.seed(2)
  n <- 20
  J <- 8
  X <- cbind(1, rep(0:1, each = n / 2))
  Y <- matrix(rpois(n * J, 5), n, J)
  Y[X[, 2] == 0, 1] <- 0
  Y[X[, 2] == 1, 2] <- 0
  closed_form <- emuFit_micro_penalized(X, Y, verbose = FALSE)
  iterative <- emuFit_micro_penalized(X, Y, verbose = FALSE, use_discrete = FALSE,
                                      tolerance = 1e-8)
  expect_true(closed_form$convergence)
  expect_lt(max(abs(closed_form$B - iterative$B)), 1e-6)
  expect_lt(max(abs(closed_form$Y_augmented - iterative$Y_augmented)), 1e-6)
  # 1/2 added to each covariate-pattern-by-category total
  expect_equal(unname(closed_form$B[2, 1] - closed_form$B[2, 3]),
               log((sum(Y[X[, 2] == 1, 1]) + 0.5) / 0.5) -
                 log((sum(Y[X[, 2] == 1, 3]) + 0.5) / (sum(Y[X[, 2] == 0, 3]) + 0.5)))
})

test_that("closed-form discrete penalized fit matches iterative fit, unbalanced three-level design with shuffled rows", {
  set.seed(3)
  J <- 10
  g <- sample(rep(1:3, c(4, 7, 11)))
  n <- length(g)
  X <- cbind(1, g == 2, g == 3)
  Y <- matrix(rpois(n * J, 3), n, J)
  Y[g == 3, 4] <- 0
  constraint_fn <- rep(list(function(x) mean(x)), 3)
  closed_form <- emuFit_micro_penalized(X, Y, constraint_fn = constraint_fn,
                                        verbose = FALSE)
  iterative <- emuFit_micro_penalized(X, Y, constraint_fn = constraint_fn,
                                      verbose = FALSE, use_discrete = FALSE,
                                      tolerance = 1e-8)
  expect_lt(max(abs(closed_form$B - iterative$B)), 1e-6)
  expect_lt(max(abs(closed_form$Y_augmented - iterative$Y_augmented)), 1e-6)
})

test_that("closed-form discrete penalized fit errors on rank-deficient design", {
  # three distinct rows, third column twice the second
  X <- cbind(1, c(0, 0, 1, 1, 2), c(0, 0, 2, 2, 4))
  Y <- matrix(1:15, 5, 3)
  expect_error(emuFit_micro_penalized(X, Y, verbose = FALSE), "rank-deficient")
})

test_that("emuFit with estimates_only returns the same estimates without test_kj", {
  set.seed(4)
  n <- 20
  J <- 12
  X <- cbind(1, rep(0:1, each = n / 2))
  colnames(X) <- c("(Intercept)", "group")
  Y <- matrix(rpois(n * J, 5), n, J)
  Y[1, colSums(Y) == 0] <- 1
  colnames(Y) <- paste0("cat", 1:J)
  est_only <- emuFit(Y = Y, X = X, estimates_only = TRUE)
  full <- emuFit(Y = Y, X = X, run_score_tests = FALSE, compute_cis = TRUE)
  expect_equal(est_only$B, full$B)
  expect_equal(est_only$z_hat, full$z_hat)
  expect_equal(est_only$coef$estimate, full$coef$estimate)
  expect_equal(est_only$coef$category, colnames(Y))
  expect_equal(est_only$coef$covariate, rep("group", J))
  expect_null(est_only$I)
  expect_true(all(is.na(est_only$coef$pval)))
  expect_warning(emuFit(Y = Y, X = X, estimates_only = TRUE, run_score_tests = TRUE),
                 "estimates_only")
})

test_that("emuFit with estimates_only matches full fit for a continuous covariate", {
  set.seed(5)
  n <- 15
  J <- 6
  X <- cbind(1, rnorm(n))
  Y <- matrix(rpois(n * J, 5), n, J)
  Y[1, colSums(Y) == 0] <- 1
  est_only <- emuFit(Y = Y, X = X, estimates_only = TRUE)
  full <- emuFit(Y = Y, X = X, run_score_tests = FALSE, compute_cis = TRUE)
  expect_equal(est_only$B, full$B)
  expect_equal(est_only$z_hat, full$z_hat)
})

test_that("emuFit with refit = FALSE uses closed-form augmentations for discrete designs", {
  set.seed(6)
  n <- 20
  J <- 8
  X <- cbind(1, rep(0:1, each = n / 2))
  Y <- matrix(rpois(n * J, 5), n, J)
  Y[1, colSums(Y) == 0] <- 1
  fit <- emuFit(Y = Y, X = X, estimates_only = TRUE)
  refit <- emuFit(Y = Y, X = X, B = fit$B, refit = FALSE, estimates_only = TRUE)
  expect_lt(max(abs(refit$Y_augmented - fit$Y_augmented)), 1e-10)
  expect_lt(max(abs(refit$z_hat - fit$z_hat)), 1e-10)
})
