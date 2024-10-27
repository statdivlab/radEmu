
test_that("wald test gives semi-reasonable output", {
  set.seed(343234)
  n <- 20
  X <- cbind(1,rep(c(0,1),each = n/2))
  J <- 10
  b0 <- rnorm(10)
  b1 <- 1:10
  b1 <- b1 - mean(b1)
  b1[5] <- pseudohuber_center(b1[-5],0.1)
  b0 <- b0 - mean(b0)
  Y <- radEmu:::simulate_data(n = n,
                              J = J,
                              X = X,
                              b0 = b0,
                              b1 = b1,
                              distn = "ZINB",
                              zinb_size = 3,
                              zinb_zero_prop = 0.6,
                              mean_z = 8)
  
  k_constr <- 2
  j_constr <- 5
  p <- 2
  
  constraint_fn <- function(x){ pseudohuber_center(x,0.1)}
  
  ##### Arguments to fix:
  
  constraint_grad_fn <- function(x){dpseudohuber_center_dx(x,0.1)}
  
  constraint_hess_fn <- function(x,ind_1,ind_2){hess_pseudohuber_center(x,0.1,ind_1,ind_2)}
  b[2,4] <- constraint_fn(b[2,-4])
  
  X_cup <- X_cup_from_X(X,J)
  
  full_fit <- emuFit_micro_penalized(X = X,
                                     Y = Y,
                                     B = NULL,
                                     constraint_fn = constraint_fn,
                                     tolerance = 1e-3,
                                     verbose = FALSE)
  
  wald_result <-
    micro_wald(Y = full_fit$Y,
               X,
               X_cup = X_cup,
               B = full_fit$B,
               test_kj = data.frame(k = 2, j = 4),
               constraint_fn = constraint_fn,
               constraint_grad_fn = constraint_grad_fn,
               nominal_coverage = 0.95)
  
  expect_true(is.data.frame(wald_result$coefficients))
  expect_true(wald_result$coefficients$pval>0.1)
  expect_true(is.list(wald_result))
  expect_true(ncol(wald_result$I) ==20)
  expect_equal(wald_result$coefficients$pval, 0.61, tolerance = 0.02)
  
})



test_that("wald test gives semi-reasonable output with continuous covariate", {
  set.seed(9944234)
  n <- 20
  X <- cbind(1,rnorm(n))
  J <- 10
  z <- rnorm(n) +8
  b0 <- rnorm(10)
  b1 <- 1:10
  b1 <- b1 - mean(b1)
  b1[5] <- pseudohuber_center(b1[-5],0.1)
  b0 <- b0 - mean(b0)
  Y <- radEmu:::simulate_data(n = n,
                              J = J,
                              X = X,
                              b0 = b0,
                              b1 = b1,
                              distn = "ZINB",
                              zinb_size = 3,
                              zinb_zero_prop = 0.6,
                              mean_z = 8)
  
  k_constr <- 2
  j_constr <- 5
  p <- 2
  
  constraint_fn <- function(x){ pseudohuber_center(x,0.1)}
  
  ##### Arguments to fix:
  
  constraint_grad_fn <- function(x){dpseudohuber_center_dx(x,0.1)}
  
  constraint_hess_fn <- function(x,ind_1,ind_2){hess_pseudohuber_center(x,0.1,ind_1,ind_2)}
  b[2,4] <- constraint_fn(b[2,-4])
  
  X_cup <- X_cup_from_X(X,J)
  
  full_fit <- emuFit_micro_penalized(X = X,
                                     Y = Y,
                                     B = NULL,
                                     constraint_fn = constraint_fn,
                                     tolerance = 1e-3,
                                     verbose = FALSE)
  
  wald_result <- micro_wald(Y = full_fit$Y,
                            X,
                            X_cup = X_cup,
                            B = full_fit$B,
                            test_kj = data.frame(k = 2, j = 4),
                            constraint_fn = constraint_fn,
                            constraint_grad_fn = constraint_grad_fn,
                            nominal_coverage = 0.95)
  
  wald_result_for_an_alternative <- micro_wald(Y = full_fit$Y,
                                               X,
                                               X_cup = X_cup,
                                               B = full_fit$B,
                                               test_kj = data.frame(k = 2, j = 10),
                                               constraint_fn = constraint_fn,
                                               constraint_grad_fn = constraint_grad_fn,
                                               nominal_coverage = 0.95)
  
  
  expect_true(is.data.frame(wald_result$coefficients))
  expect_true(is.list(wald_result))
  expect_true(ncol(wald_result$I) ==20)
  expect_equal(wald_result$coefficients$pval, 0.11, tolerance = 0.03)
  
  expect_true(wald_result_for_an_alternative$coefficients$pval < 0.01)
  
})
