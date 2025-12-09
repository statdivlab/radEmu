
test_that("new discrete is correct", {
  
  # set.seed(1)
  n <- 40
  J <- 30
  p <- 5
  beta <- matrix(rnorm(p*J, mean = 2, sd=1), ncol = J)
  X <- cbind(1, 
             c(rep(0, n/2), rep(1, n/2)), 
             c(rep(0, 5*n/8), rep(1, 3*n/8)), 
             c(rep(0, 3*n/4), rep(1, n/4)),
             c(rep(0, 7*n/8), rep(1, n/8)))
  Y <- matrix(rpois(n=n*J, lambda=exp(1 + X %*% beta)), ncol = J)
  Y <- pmax(Y, 1)
  
  my_jstar <- 9
  my_kstar <- 4
  my_jref <- J
  
  t_discrete <- system.time({
    out_discrete <- fit_null_discrete_pseudohuber(Y = Y, X = X, k_constr = my_kstar, j_constr = my_jstar, j_ref = my_jref, 
                                                  tol = 1e-12)
  })
  
  ## check that the (my_kstar, my_kstar) constraint is satisfied
  expect_equal(pseudohuber_median(out_discrete$B[my_kstar, -my_jstar]), 
               out_discrete$B[my_kstar, my_jstar] )
  
  ## check that the j_ref constraint is satisfied
  expect_equal(out_discrete$B[, my_jref], 
               rep(0, p))
  
  ## check permutation invariant 
  X2 <- X[, c(1, (p:2))]
  new_order <- c(sample(1:(J-1), J-1), J)
  new_j_star <- which(new_order == my_jstar)
  Y2 <- Y[, new_order]
  t_discrete2 <- system.time({
    out_discrete2 <- fit_null_discrete_pseudohuber(Y = Y2, X = X2, k_constr = which(c(1, (p:2)) == my_kstar), 
                                                   j_constr = new_j_star, j_ref = my_jref, 
                                                   tol = 1e-12)
  })
  
  z_discrete <- update_z(Y, X, out_discrete$B)
  log_means_discrete <- X %*% out_discrete$B + matrix(z_discrete, ncol = 1) %*% matrix(1, ncol = J, nrow = 1)
  ll_new_discrete <- sum(Y * log_means_discrete - exp(log_means_discrete))
  
  z_discrete2 <- update_z(Y2, X2, out_discrete2$B)
  log_means_discrete2 <- X2 %*% out_discrete2$B + matrix(z_discrete2, ncol = 1) %*% matrix(1, ncol = J, nrow = 1)
  ll_new_discrete2 <- sum(Y2 * log_means_discrete2 - exp(log_means_discrete2))
  
  ## expect them to be close in estimates and likelihood
  expect_true(max(abs(out_discrete$B[, new_order] - out_discrete2$B[c(1, (p:2)), ])) < 1e-4)
  expect_equal(ll_new_discrete2, ll_new_discrete)
  
})



test_that("new discrete is correct with flexible jref and jstar", {
  
  # set.seed(3)
  n <- 24
  J <- 10
  p <- 5
  beta <- matrix(rnorm(p*J, mean = 1, sd=1), ncol = J)
  X <- cbind(1, 
             c(rep(0, n/2), rep(1, n/2)), 
             c(rep(0, 5*n/8), rep(1, 3*n/8)), 
             c(rep(0, 3*n/4), rep(1, n/4)),
             c(rep(0, 7*n/8), rep(1, n/8)))
  Y <- matrix(rpois(n=n*J, lambda=exp(1 + X %*% beta)), ncol = J)
  Y <- pmax(Y, 1)
  
  my_jstar <- 9
  my_jref <- 1
  my_kstar <- 4
  
  out_discrete91 <- fit_null_discrete_pseudohuber(Y = Y, X = X, k_constr = my_kstar, j_constr = my_jstar, j_ref = my_jref, 
                                                  tol = 1e-12)
  
  ## check that the (my_kstar, my_kstar) constraint is satisfied
  expect_equal(pseudohuber_median(out_discrete91$B[my_kstar, -my_jstar]), 
               out_discrete91$B[my_kstar, my_jstar] )
  
  ## check that the j_ref constraint is satisfied
  expect_equal(out_discrete91$B[, my_jref], 
               rep(0, p))
  
  my_jstar <- 4
  my_jref <- 3
  
  out_discrete43 <- fit_null_discrete_pseudohuber(Y = Y, X = X, k_constr = my_kstar, j_constr = my_jstar, j_ref = my_jref, 
                                                  tol = 1e-12)
  
  ## check that the (my_kstar, my_kstar) constraint is satisfied
  expect_equal(pseudohuber_median(out_discrete43$B[my_kstar, -my_jstar]), 
               out_discrete43$B[my_kstar, my_jstar] )
  
  ## check that the j_ref constraint is satisfied
  expect_equal(out_discrete43$B[, my_jref], 
               rep(0, p))
  
  my_jstar <- 4
  my_jref <- 6
  
  out_discrete46 <- fit_null_discrete_pseudohuber(Y = Y, X = X, k_constr = my_kstar, j_constr = my_jstar, j_ref = my_jref, 
                                                  tol = 1e-12)
  
  ## should be same value across rows
  expect_equal((out_discrete43$B - out_discrete46$B) %>% round(4) %>% apply(1, sd) %>% max, 0)
  
})




test_that("new discrete is fast", {
  
  skip("Too slow for automated testing; 20 mins on laptop")
  
  set.seed(2)
  n <- 12
  J <- 10
  p <- 3
  beta <- matrix(rnorm(p*J, mean = 2, sd=1), ncol = J)
  X <- cbind(1, c(rep(0, n/2), rep(1, n/2)), c(rep(0, 3*n/4), rep(1, n/4)))
  Y <- matrix(rpois(n=n*J, lambda=exp(2 + X %*% beta)), ncol = J)
  Y <- pmax(Y, 1)
  
  my_jstar <- J-1
  my_kstar <- 2
  my_jref <- J
  
  t_discrete <- system.time({
    out_discrete <- fit_null_discrete_pseudohuber(Y = Y, X = X, k_constr = my_kstar, j_constr = my_jstar, j_ref = my_jref, 
                                                  tol = 1e-12)
  })
  
  t_orig <- system.time({
    out_orig <- fit_null(B=matrix(0, nrow = p, ncol = J), 
                         Y=Y, X = X, 
                         k_constr=my_kstar, j_constr=my_jstar, j_ref=my_jref,
                         constraint_fn=list(pseudohuber_median, pseudohuber_median),
                         constraint_grad_fn=list(radEmu::dpseudohuber_median_dx, radEmu::dpseudohuber_median_dx),
                         B_tol=1e-8, constraint_tol=1e-10, maxit = 1e8)
  }) #  816.931  489.308 1317.577 
  
  ## check faster
  expect_lt(t_discrete[3], t_orig[3])
  
  ## check that the (my_kstar, my_kstar) constraint is satisfied
  expect_equal(pseudohuber_median(out_discrete$B[my_kstar, -my_jstar]), 
               out_discrete$B[my_kstar, my_jstar] )
  
  expect_equal(pseudohuber_median(out_orig$B[my_kstar, -my_jstar]), 
               out_orig$B[my_kstar, my_jstar] )
  
  expect_lte(max(abs(out_discrete$B - out_orig$B)), 1e-3)
  
  z_discrete <- update_z(Y, X, out_discrete$B)
  log_means_discrete <- X %*% out_discrete$B + matrix(z_discrete, ncol = 1) %*% matrix(1, ncol = J, nrow = 1)
  ll_new_discrete <- sum(Y * log_means_discrete - exp(log_means_discrete))
  
  z_orig <- update_z(Y, X, out_orig$B)
  log_means_orig <- X %*% out_orig$B + matrix(z_orig, ncol = 1) %*% matrix(1, ncol = J, nrow = 1)
  ll_new_orig <- sum(Y * log_means_orig - exp(log_means_orig))
  
  ## expect the new likelihood is no worse than 1% lower
  expect_true(ll_new_discrete > ll_new_orig*0.99)
  
})


test_that("new discrete null aligns with older code", {
  
  skip("Too slow for automated testing; 3 mins on laptop")
  
  set.seed(4)
  n <- 12
  J <- 5
  p <- 3
  beta <- matrix(rnorm(p*J, mean = 2, sd=1), ncol = J)
  X <- cbind(1, c(rep(0, n/3), rep(1, n/3), rep(0, n/3)), c(rep(0, n/3), rep(0, n/3), rep(1, n/3)))
  Y <- matrix(rpois(n=n*J, lambda=exp(2 + X %*% beta)), ncol = J)
  Y <- pmax(Y, 1)
  
  my_jstar <- J-2
  my_jref <- J-1
  my_kstar <- 2
  
  t_discrete <- system.time({
    out_discrete <- fit_null_discrete_pseudohuber(Y = Y, X = X, k_constr = my_kstar, j_constr = my_jstar, j_ref = my_jref, 
                                                  tol = 1e-12)
  })
  
  t_orig <- system.time({
    out_orig <- fit_null(B=matrix(0, nrow = p, ncol = J), 
                         Y=Y, X = X, 
                         k_constr=my_kstar, j_constr=my_jstar, j_ref=my_jref,
                         constraint_fn=list(pseudohuber_median, pseudohuber_median),
                         constraint_grad_fn=list(radEmu::dpseudohuber_median_dx, radEmu::dpseudohuber_median_dx),
                         B_tol=1e-3, constraint_tol=1e-6, maxit = 1e8)
  }) # 94.249  66.621 163.194 
  
  ## check faster
  expect_lt(t_discrete[3], t_orig[3])
  
  ## check align
  expect_lte(max(abs(out_discrete$B - out_orig$B)), 0.02)
  
  ## check no worse
  z_discrete <- update_z(Y, X, out_discrete$B)
  log_means_discrete <- X %*% out_discrete$B + matrix(z_discrete, ncol = 1) %*% matrix(1, ncol = J, nrow = 1)
  ll_new_discrete <- sum(Y * log_means_discrete - exp(log_means_discrete))
  
  z_orig <- update_z(Y, X, out_orig$B)
  log_means_orig <- X %*% out_orig$B + matrix(z_orig, ncol = 1) %*% matrix(1, ncol = J, nrow = 1)
  ll_new_orig <- sum(Y * log_means_orig - exp(log_means_orig))
  
  ## expect the new likelihood is no worse than 1% lower
  expect_true(ll_new_discrete > ll_new_orig*0.99)
  
})