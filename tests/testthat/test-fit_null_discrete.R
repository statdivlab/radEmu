
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
    out_discrete <- fit_null_discrete(Y = Y, X = X, k_constr = my_kstar, j_constr = my_jstar, 
                                      j_ref = my_jref, tol = 1e-12)
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
    out_discrete2 <- fit_null_discrete(Y = Y2, X = X2, k_constr = which(c(1, (p:2)) == my_kstar), 
                                       j_constr = new_j_star, j_ref = my_jref, tol = 1e-12)
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


# tests to compare to fit_null
test_that("new discrete is correct with flexible jref and jstar", {
  
  # set.seed(3)
  n <- 80 # increasing n here so that this should hold across random variation
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
  
  out_discrete91 <- fit_null_discrete(Y = Y, X = X, k_constr = my_kstar, j_constr = my_jstar, j_ref = my_jref, 
                                      tol = 1e-12)
  
  ## check that the (my_kstar, my_kstar) constraint is satisfied
  expect_equal(pseudohuber_median(out_discrete91$B[my_kstar, -my_jstar]), 
               out_discrete91$B[my_kstar, my_jstar] )
  
  ## check that the j_ref constraint is satisfied
  expect_equal(out_discrete91$B[, my_jref], 
               rep(0, p))
  
  my_jstar <- 4
  my_jref <- 3
  
  out_discrete43 <- fit_null_discrete(Y = Y, X = X, k_constr = my_kstar, j_constr = my_jstar, j_ref = my_jref, 
                                      tol = 1e-12)
  
  ## check that the (my_kstar, my_kstar) constraint is satisfied
  expect_equal(pseudohuber_median(out_discrete43$B[my_kstar, -my_jstar]), 
               out_discrete43$B[my_kstar, my_jstar] )
  
  ## check that the j_ref constraint is satisfied
  expect_equal(out_discrete43$B[, my_jref], 
               rep(0, p))
  
  my_jstar <- 4
  my_jref <- 6
  
  out_discrete46 <- fit_null_discrete(Y = Y, X = X, k_constr = my_kstar, j_constr = my_jstar, j_ref = my_jref, 
                                      tol = 1e-12)
  
  ## should be same value across rows
  # commenting this test out because it is typically true but not always across all
  # random draws
  # expect_equal((out_discrete43$B - out_discrete46$B) %>% round(4) %>% apply(1, sd) %>% max, 0)
  
})

test_that("new discrete is fast", {
  
  #skip("Too slow for automated testing; 20 mins on laptop")
  
  set.seed(2)
  n <- 60
  J <- 10
  p <- 3
  beta <- matrix(rnorm(p*J, mean = 2, sd=1), ncol = J)
  X <- cbind(1, c(rep(0, n/2), rep(1, n/2)), c(rep(0, 3*n/4), rep(1, n/4)))
  Y <- matrix(rpois(n=n*J, lambda=exp(2 + X %*% beta)), ncol = J)
  Y <- pmax(Y, 1)
  
  my_jstar <- J-1
  my_kstar <- 2
  my_jref <- J
  
  alt_fit <- emuFit(Y = Y, X = X, run_score_tests = FALSE, compute_cis = FALSE,
                    match_row_names = FALSE)
  
  t_discrete <- system.time({
    out_discrete <- fit_null_discrete(Y = Y, X = X, k_constr = my_kstar, j_constr = my_jstar, j_ref = my_jref)
  })
  
  t_orig <- system.time({
    out_orig <- fit_null(B=alt_fit$B, 
                         Y=Y, X = X, 
                         k_constr=my_kstar, j_constr=my_jstar, j_ref=my_jref,
                         constraint_fn=list(pseudohuber_median, pseudohuber_median),
                         constraint_grad_fn=list(radEmu::dpseudohuber_median_dx, radEmu::dpseudohuber_median_dx))
  }) 
  
  ## check faster
  expect_lt(t_discrete[3], t_orig[3])
  
  ## check that the (my_kstar, my_kstar) constraint is satisfied
  expect_equal(pseudohuber_median(out_discrete$B[my_kstar, -my_jstar]), 
               out_discrete$B[my_kstar, my_jstar] )
  
  expect_equal(pseudohuber_median(out_orig$B[my_kstar, -my_jstar]), 
               out_orig$B[my_kstar, my_jstar] , tolerance = 1e-5)
  
  #expect_lte(max(abs(out_discrete$B - out_orig$B)), 1e-3)
  
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
  
  skip("Too slow for automated testing; ~3 mins on laptop")
  
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
    out_discrete <- fit_null_discrete(Y = Y, X = X, k_constr = my_kstar, j_constr = my_jstar, j_ref = my_jref, 
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

# tests to compare to fit_null_symmetric
set.seed(59542234)
n <- 40
J <- 20
X <- cbind(1, rep(c(0, 1), each = n / 2))
b0 <- rnorm(J)
b1 <- seq(1, 10, length.out = J)
b1 <- b1 - mean(b1)
b0 <- b0 - mean(b0)
Y <- radEmu:::simulate_data(
  n = n,
  J = J,
  X = X,
  b0 = b0,
  b1 = b1,
  distn = "Poisson",
  mean_z = 8
)

test_that("we can choose null fit algorithm argument", {
  fit1 <- emuFit(X = X, Y = Y, test_kj = data.frame(k = 2, j = 5), null_fit_alg = "constraint_sandwich", tolerance = 1e-2,
                 match_row_names = FALSE, null_diagnostic_plots = TRUE)
  fit2 <- emuFit(X = X, Y = Y, test_kj = data.frame(k = 2, j = 5), null_fit_alg = "augmented_lagrangian", tolerance = 1e-2,
                 match_row_names = FALSE, null_diagnostic_plots = TRUE)
  fit3 <- emuFit(X = X, Y = Y, test_kj = data.frame(k = 2, j = 5), null_fit_alg = "discrete", tolerance = 1e-2,
                 match_row_names = FALSE, null_diagnostic_plots = TRUE)
  
  # expect pretty similar results
  expect_true(all.equal(fit1$coef$score_stat[5], fit2$coef$score_stat[5], tol = 0.01))
  expect_true(all.equal(fit1$coef$score_stat[5], fit3$coef$score_stat[5], tol = 0.01))
  # but not exactly the same, because different algorithms are being used
  expect_false(fit1$coef$score_stat[5] == fit3$coef$score_stat[5])
})

test_that("compare timing", {
  
  skip("long for automatic testing")
  
  sand_start <- proc.time()
  fit1 <- emuFit(X = X, Y = Y, test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "constraint_sandwich", tolerance = 1e-2,
                 match_row_names = FALSE, null_diagnostic_plots = T)
  sand_end <- proc.time() - sand_start
  discrete_start <- proc.time()
  fit2 <- emuFit(X = X, Y = Y, test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "discrete", tolerance = 1e-2,
                 match_row_names = FALSE, null_diagnostic_plots = T, maxit_null = 3000)
  discrete_end <- proc.time() - discrete_start
  lag_start <- proc.time()
  fit3 <- emuFit(X = X, Y = Y, test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "augmented_lagrangian", tolerance = 1e-2,
                 match_row_names = FALSE, null_diagnostic_plots = T, maxit_null = 3000)
  lag_end <- proc.time() - lag_start
  
  sand_end; discrete_end; lag_end
  lik1 <- sapply(fit1$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
  lik2 <- sapply(fit2$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
  lik3 <- sapply(fit3$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
  summary(lik1 - lik2)
  summary(lik3 - lik2)
  summary(lik1 - lik3)
  plot(lik1 - lik2, fit1$coef$score_stat - fit2$coef$score_stat)
  plot(lik3 - lik2, fit3$coef$score_stat - fit2$coef$score_stat)
  plot(lik1 - lik3, fit1$coef$score_stat - fit3$coef$score_stat)
  
  # sandwich is 8 seconds, discrete 11, augmented lagrangian 52
  # when score is different, discrete has higher ll
  # ll varies between methods but more often discrete has higher ll 
})

test_that("compare timing, ZINB", {
  
  skip("long for automatic testing")
  
  Y_zinb <- radEmu:::simulate_data(
    n = n,
    J = J,
    X = X,
    b0 = b0,
    b1 = b1,
    distn = "ZINB",
    mean_z = 8,
    zinb_size = 5,
    zinb_zero_prop = 0.6
  )
  
  sand_start <- proc.time()
  fit1 <- emuFit(X = X, Y = Y_zinb, test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "constraint_sandwich", tolerance = 1e-2,
                 match_row_names = FALSE, null_diagnostic_plots = T, verbose = TRUE)
  sand_end <- proc.time() - sand_start
  discrete_start <- proc.time()
  fit2 <- emuFit(X = X, Y = Y_zinb, test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "discrete", tolerance = 1e-2,
                 match_row_names = FALSE, null_diagnostic_plots = T, 
                 maxit_null = 4000, verbose = TRUE)
  discrete_end <- proc.time() - discrete_start
  lag_start <- proc.time()
  fit3 <- emuFit(X = X, Y = Y_zinb, test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "augmented_lagrangian", tolerance = 1e-2,
                 match_row_names = FALSE, null_diagnostic_plots = T, verbose = TRUE)
  lag_end <- proc.time() - lag_start
  
  lik1 <- sapply(fit1$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
  lik2 <- sapply(fit2$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
  lik3 <- sapply(fit3$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
  summary(lik1 - lik2)
  summary(lik3 - lik2)
  summary(lik1 - lik3)
  plot(lik1 - lik2, fit1$coef$score_stat - fit2$coef$score_stat)
  plot(lik3 - lik2, fit3$coef$score_stat - fit2$coef$score_stat)
  plot(lik1 - lik3, fit1$coef$score_stat - fit3$coef$score_stat)
  
  sand_end; discrete_end; lag_end
  # discrete is 6 seconds, sandwich is 12 seconds, augmented lagrangian in 50 seconds
  # ll is nearly always larger for discrete than other algorithms
})

# the next set of tests compare the timing of fit_null_discrete to fit_null_symmetric
# for a variety of n, J, and p, using the soil dataset included in `corncob` and the
# wirbel dataset included in `radEmu`. Each example runs either 10, 20, or 30 robust 
# score tests and compares across the two methods

# tldr:

# wirbel
# n = 126, J = 128, p = 2, sandwich 42 seconds, discrete 16 seconds
# n = 566, J = 133, p = 5, sandwich 306 seconds, discrete 334 seconds
# n = 126, J = 430, p = 2, sandwich 143 seconds, discrete 421 seconds
# n = 126, J = 758, p = 2, sandwich 8 minutes, discrete 60 minutes

# soil
# n = 119, J = 109, p = 3, sandwich 35 seconds, discrete 12 seconds
# n = 119, J = 147, p = 3, sandwich 50 seconds, discrete 49 seconds
# n = 64, J = 234, p = 2, sandwich 121 seconds, discrete 165 seconds
# n = 64, J = 242, p = 2, sandwich 101 seconds, discrete 95 seconds
# n = 64, J = 479, p = 2, sandwich 140 seconds, discrete 534 seconds

# generally discrete attains a higher ll, but the ll's are quite similar in most
# cases, as are the resulting test statistics 

# overall heuristic from this investigation is to use the discrete algorithm
# when J < 150, and the sandwich constraint algorithm when J >= 150

# tests producing these results are included below

# # commenting out to avoid note about corncob data not being part of package
# # however keeping because these are useful for testing
# test_that("test corncob data, n = 64, J = 234, p = 2", {
# 
#   skip("don't want to test automatically")
# 
#   # corncob data
# 
#   library(corncob)
#   library(phyloseq)
#   data(soil_phylo_sample)
#   data(soil_phylo_otu)
#   soil_phylo <- phyloseq::phyloseq(phyloseq::sample_data(soil_phylo_sample),
#                                    phyloseq::otu_table(soil_phylo_otu, taxa_are_rows = TRUE),
#                                    phyloseq::tax_table(soil_phylo_taxa))
#   soil_phylo <- subset_samples(soil_phylo, Day %in% 0:1 & Amdmt %in% 0:1)
#   soil_phylo <- tax_glom(soil_phylo, taxrank = "Genus")
#   soil_samp <- data.frame(sample_data(soil_phylo))
#   soil_otu <- t(unclass(otu_table(soil_phylo)))
# 
#   to_rm <- which(colSums(soil_otu) == 0)
#   soil_otu <- soil_otu[, -to_rm]
# 
#   soil_fit <- emuFit(formula = ~ Amdmt,
#                      data = soil_samp,
#                      Y = soil_otu,
#                      run_score_tests = FALSE, compute_cis = FALSE,
#                      tolerance = 1e-5, verbose = "development")
#   start_sand <- proc.time()
#   corn_sand1_new <- emuFit(formula = ~ Amdmt,
#                        data = soil_samp,
#                        Y = soil_otu,
#                        fitted_model = soil_fit,
#                        refit = FALSE,
#                        compute_cis = FALSE,
#                        test_kj = data.frame(k = 2, j = 1:30), null_fit_alg = "constraint_sandwich",
#                        verbose = TRUE, null_diagnostic_plots = T)
#   end_sand <- proc.time() - start_sand
#   # 3, 5, 3, 1, 4, 3, 4, 4, 8, 3, 3, 8, 5, 3, 3, 3, 2, 3, 2, 1, 4, 4, 3, 1, 6, 2, 4, 4, 5
#   start_discrete <- proc.time()
#   corn_discrete1 <- emuFit(formula = ~ Amdmt,
#                       data = soil_samp,
#                       Y = soil_otu,
#                       fitted_model = soil_fit,
#                       refit = FALSE,
#                       compute_cis = FALSE,
#                       test_kj = data.frame(k = 2, j = 1:30), null_fit_alg = "discrete",
#                       verbose = TRUE, null_diagnostic_plots = T)
#   end_discrete <- proc.time() - start_discrete
#   # 4, 13, 5, 9, 10, 9, 12, 10, 16, 7, 13, 34, 7, 11, 9, 5, 3, 9, 9, 5, 8, 8, 11, 14, 5, 22, 6, 12, 10, 8
#   
#   end_sand[3]; end_discrete[3]
# 
#   lik_sand <- sapply(corn_sand1_new$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
#   lik_discrete <- sapply(corn_discrete1$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
#   plot(lik_discrete - lik_sand, corn_sand1_new$coef$score_stat[1:30] - corn_discrete1$coef$score_stat[1:30])
# 
#   summary(lik_discrete - lik_sand)
#   summary(corn_sand1_new$coef$score_stat[1:30] - corn_discrete1$coef$score_stat[1:30])
# 
#   # n = 64, J = 234, p = 2
#   # sandwich takes 121 seconds, discrete takes 165 seconds
#   # ll is nearly always higher (sometimes significantly) for discrete
# })

test_that("test wirbel data, n = 126, J = 128, p = 2", {
  
  skip("don't want to test automatically")
  
  data("wirbel_sample")
  data("wirbel_otu")
  data("wirbel_taxonomy")
  wirbel_phylo <- phyloseq::phyloseq(phyloseq::sample_data(wirbel_sample),
                                     phyloseq::otu_table(wirbel_otu, taxa_are_rows = FALSE),
                                     phyloseq::tax_table(wirbel_taxonomy))
  wirbel_genus <- phyloseq::tax_glom(wirbel_phylo, taxrank = "genus")
  wirbel_genus_ch <- phyloseq::subset_samples(wirbel_genus, Country == "CHI")
  zero_taxa <- phyloseq::taxa_sums(wirbel_genus_ch) == 0
  wirbel_genus_ch <- phyloseq::prune_taxa(zero_taxa == FALSE, wirbel_genus_ch)
  ch_fit_full <- emuFit(formula = ~ Group, 
                        Y = wirbel_genus_ch, run_score_tests = FALSE, 
                        compute_cis = FALSE, verbose = "development",
                        tolerance = 1e-5)
  
  start_sand <- proc.time()
  wirb_sand <- emuFit(formula = ~ Group, 
                      Y = wirbel_genus_ch,
                      fitted_model = ch_fit_full,
                      refit = FALSE, 
                      compute_cis = FALSE, 
                      test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "constraint_sandwich",
                      verbose = TRUE, null_diagnostic_plots = T)
  end_sand <- proc.time() - start_sand
  # 5, 1, 1, 1, 3, 1, 3, 1, 3, 2, 1, 3, 3, 1, 1, 3, 1, 1, 1, 4
  start_discrete <- proc.time()
  wirb_discrete <- emuFit(formula = ~ Group, 
                      Y = wirbel_genus_ch,
                      fitted_model = ch_fit_full,
                      refit = FALSE, 
                      compute_cis = FALSE, 
                      test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "discrete",
                      verbose = TRUE, null_diagnostic_plots = T)
  end_discrete <- proc.time() - start_discrete
  # 3, 1, 1, 1, 2, 1, 2, 1, 2, 2, 0, 2, 2, 1, 0, 4, 1, 1, 1, 3
  end_sand[3]; end_discrete[3]
  
  lik_sand <- sapply(wirb_sand$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
  lik_discrete <- sapply(wirb_discrete$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
  plot(lik_discrete - lik_sand, wirb_sand$coef$score_stat[1:20] - wirb_discrete$coef$score_stat[1:20])
  
  summary(lik_discrete - lik_sand)
  summary(abs(wirb_sand$coef$score_stat[1:20] - wirb_discrete$coef$score_stat[1:20]))
  
  # n = 126, J = 128, p = 2
  # results: sandwich is 42 seconds, discrete is 16 seconds
  # ll is almost always higher for discrete but generally ll and test stats are similar
  
})

test_that("test wirbel data country, n = 566, J = 133, p = 5", {
  
  skip("don't want to test automatically")
  
  data("wirbel_sample")
  data("wirbel_otu")
  data("wirbel_taxonomy")
  wirbel_phylo <- phyloseq::phyloseq(phyloseq::sample_data(wirbel_sample),
                                     phyloseq::otu_table(wirbel_otu, taxa_are_rows = FALSE),
                                     phyloseq::tax_table(wirbel_taxonomy))
  wirbel_genus <- phyloseq::tax_glom(wirbel_phylo, taxrank = "genus")
  #wirbel_genus_ch <- phyloseq::subset_samples(wirbel_genus, Country == "CHI")
  zero_taxa <- phyloseq::taxa_sums(wirbel_genus) == 0
  wirbel_genus <- phyloseq::prune_taxa(zero_taxa == FALSE, wirbel_genus)
  ch_fit_full <- emuFit(formula = ~ Country, 
                        Y = wirbel_genus, run_score_tests = FALSE, 
                        compute_cis = FALSE, verbose = "development",
                        tolerance = 1e-5)
  
  start_sand <- proc.time()
  wirb_sand1 <- emuFit(formula = ~ Country, 
                      Y = wirbel_genus,
                      fitted_model = ch_fit_full,
                      refit = FALSE, 
                      compute_cis = FALSE, 
                      test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "constraint_sandwich",
                      verbose = TRUE, null_diagnostic_plots = T)
  end_sand <- proc.time() - start_sand
  # 28, 11, 11, 7, 13, 14, 17, 7, 7, 10, 7, 54, 25, 6, 11, 32, 7, 9, 13, 10
  start_discrete <- proc.time()
  wirb_discrete1 <- emuFit(formula = ~ Country, 
                          Y = wirbel_genus,
                          fitted_model = ch_fit_full,
                          refit = FALSE, 
                          compute_cis = FALSE, 
                          test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "discrete",
                          verbose = TRUE, null_diagnostic_plots = T)
  end_discrete <- proc.time() - start_discrete
  # 11, 9, 9, 6, 13, 12, 18, 6, 8, 8, 7, 27, 31, 6, 9, 10, 6, 6, 9, 16
  end_sand[3]; end_discrete[3]
  
  lik_sand <- sapply(wirb_sand1$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
  lik_discrete <- sapply(wirb_discrete1$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
  plot(lik_discrete - lik_sand, wirb_sand1$coef$score_stat[1:20] - wirb_discrete1$coef$score_stat[1:20])
  
  summary(lik_discrete - lik_sand)
  summary(abs(wirb_sand1$coef$score_stat[1:20] - wirb_discrete1$coef$score_stat[1:20]))
  
  # n = 566, J = 133, p = 5
  # results: sandwich is 306 seconds, discrete is 334
  # discrete almost always has higher ll (sometimes by a lot)
})

test_that("test wirbel data, n = 126, J = 758, p = 2", {
  
  skip("don't want to test automatically")
  
  data("wirbel_sample")
  data("wirbel_otu")
  data("wirbel_taxonomy")
  wirbel_phylo <- phyloseq::phyloseq(phyloseq::sample_data(wirbel_sample),
                                     phyloseq::otu_table(wirbel_otu, taxa_are_rows = FALSE),
                                     phyloseq::tax_table(wirbel_taxonomy))
  wirbel_ch <- phyloseq::subset_samples(wirbel_phylo, Country == "CHI")
  zero_taxa <- phyloseq::taxa_sums(wirbel_ch) == 0
  wirbel_ch <- phyloseq::prune_taxa(zero_taxa == FALSE, wirbel_ch)
  ch_fit_all <- emuFit(formula = ~ Group, 
                        Y = wirbel_ch, run_score_tests = FALSE, 
                        compute_cis = FALSE, verbose = "development",
                        tolerance = 1e-5)
  
  start_sand <- proc.time()
  wirb_sand_all <- emuFit(formula = ~ Group, 
                      Y = wirbel_ch,
                      fitted_model = ch_fit_all,
                      refit = FALSE, 
                      compute_cis = FALSE, 
                      test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "constraint_sandwich",
                      verbose = TRUE, null_diagnostic_plots = T)
  end_sand <- proc.time() - start_sand
  # 16, 66, 20, 17, 19, 17, 17, 18, 28, 24, 19, 46, 16, 47, 17, 18, 17, 33, 19, 19
  start_discrete <- proc.time()
  wirb_discrete <- emuFit(formula = ~ Group, 
                          Y = wirbel_ch,
                          fitted_model = ch_fit_all,
                          refit = FALSE, 
                          compute_cis = FALSE, 
                          test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "discrete",
                          verbose = TRUE, null_diagnostic_plots = T)
  end_discrete <- proc.time() - start_discrete
  # 132, 11 minutes, 129, 133, 129, 160, 129, 128, 222, 139, 142, 170, 130, 6 minutes, 155, 131, 143, 211, 131, 131
  end_sand[3]; end_discrete[3]
  
  lik_sand <- sapply(wirb_sand_all$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
  lik_discrete <- sapply(wirb_discrete$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
  plot(lik_discrete - lik_sand, wirb_sand_all$coef$score_stat[1:20] - wirb_discrete$coef$score_stat[1:20])
  
  summary(lik_discrete - lik_sand)
  summary(abs(wirb_sand_all$coef$score_stat[1:20] - wirb_discrete$coef$score_stat[1:20]))
  
  # n = 126, J = 758, p = 2
  # results: sandwich is 8 minutes, discrete is 60 minutes
  # ll is almost always higher for discrete but generally ll and test stats are similar
  
})

test_that("test wirbel data, n = 126, J = 430, p = 2", {
  
  skip("don't want to test automatically")
  
  data("wirbel_sample")
  data("wirbel_otu")
  data("wirbel_taxonomy")
  wirbel_phylo <- phyloseq::phyloseq(phyloseq::sample_data(wirbel_sample),
                                     phyloseq::otu_table(wirbel_otu, taxa_are_rows = FALSE),
                                     phyloseq::tax_table(wirbel_taxonomy))
  wirbel_ch <- phyloseq::subset_samples(wirbel_phylo, Country == "CHI")
  zero_taxa <- phyloseq::taxa_sums(wirbel_ch) == 0
  wirbel_ch <- phyloseq::prune_taxa(zero_taxa == FALSE, wirbel_ch)
  wirbel_ch_clost <- phyloseq::subset_taxa(wirbel_ch, class == "Clostridia")
  ch_fit_clost <- emuFit(formula = ~ Group, 
                       Y = wirbel_ch_clost, run_score_tests = FALSE, 
                       compute_cis = FALSE, verbose = "development",
                       tolerance = 1e-5)
  
  start_sand <- proc.time()
  wirb_sand_all <- emuFit(formula = ~ Group, 
                          Y = wirbel_ch_clost,
                          fitted_model = ch_fit_clost,
                          refit = FALSE, 
                          compute_cis = FALSE, 
                          test_kj = data.frame(k = 2, j = 1:10), null_fit_alg = "constraint_sandwich",
                          verbose = TRUE, null_diagnostic_plots = T)
  end_sand <- proc.time() - start_sand
  # 7, 30, 11, 12, 21, 11, 8, 13, 17, 13
  start_discrete <- proc.time()
  wirb_discrete <- emuFit(formula = ~ Group, 
                          Y = wirbel_ch_clost,
                          fitted_model = ch_fit_clost,
                          refit = FALSE, 
                          compute_cis = FALSE, 
                          test_kj = data.frame(k = 2, j = 1:10), null_fit_alg = "discrete",
                          verbose = TRUE, null_diagnostic_plots = T)
  end_discrete <- proc.time() - start_discrete
  # 30, 92, 46, 27, 81, 33, 26, 28, 26, 33
  end_sand[3]; end_discrete[3]
  # sandwich takes 143 seconds, discrete takes 421
  
  lik_sand <- sapply(wirb_sand_all$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
  lik_discrete <- sapply(wirb_discrete$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
  plot(lik_discrete - lik_sand, wirb_sand_all$coef$score_stat[1:10] - wirb_discrete$coef$score_stat[1:10])
  
  summary(lik_discrete - lik_sand)
  summary(abs(wirb_sand_all$coef$score_stat[1:10] - wirb_discrete$coef$score_stat[1:10]))
  
  # n = 126, J = 430, p = 2
  # results: sandwich is 143 seconds, discrete is 421 seconds
  # ll is almost always higher for discrete (up to 75 higher) but generally ll and test stats are similar
  
})

# # commenting out to avoid note about corncob data not being part of package
# # however keeping because these are useful for testing
# test_that("test corncob data, n = 64, J = 479, p = 2", {
# 
#   skip("don't want to test automatically")
# 
#   # corncob data
# 
#   library(corncob)
#   library(phyloseq)
#   data(soil_phylo_sample)
#   data(soil_phylo_otu)
#   soil_phylo <- phyloseq::phyloseq(phyloseq::sample_data(soil_phylo_sample),
#                                    phyloseq::otu_table(soil_phylo_otu, taxa_are_rows = TRUE),
#                                    phyloseq::tax_table(soil_phylo_taxa))
#   soil_phylo <- subset_samples(soil_phylo, Day %in% 0:1 & Amdmt %in% 0:1)
#   soil_phylo_bact <- subset_taxa(soil_phylo, Phylum == "Bacteroidetes")
#   soil_samp <- data.frame(sample_data(soil_phylo_bact))
#   soil_otu <- t(unclass(otu_table(soil_phylo_bact)))
# 
#   to_rm <- which(colSums(soil_otu) == 0)
#   soil_otu <- soil_otu[, -to_rm]
# 
#   soil_fit <- emuFit(formula = ~ Amdmt,
#                      data = soil_samp,
#                      Y = soil_otu,
#                      run_score_tests = FALSE, compute_cis = FALSE,
#                      tolerance = 1e-5, verbose = "development")
#   start_sand <- proc.time()
#   corn_sand1_new <- emuFit(formula = ~ Amdmt,
#                        data = soil_samp,
#                        Y = soil_otu,
#                        fitted_model = soil_fit,
#                        refit = FALSE,
#                        compute_cis = FALSE,
#                        test_kj = data.frame(k = 2, j = 1:10), null_fit_alg = "constraint_sandwich",
#                        verbose = TRUE, null_diagnostic_plots = T)
#   end_sand <- proc.time() - start_sand
#   # 13, 15, 12, 6, 36, 6, 8, 18, 12, 13
#   start_discrete <- proc.time()
#   corn_discrete1 <- emuFit(formula = ~ Amdmt,
#                       data = soil_samp,
#                       Y = soil_otu,
#                       fitted_model = soil_fit,
#                       refit = FALSE,
#                       compute_cis = FALSE,
#                       test_kj = data.frame(k = 2, j = 1:10), null_fit_alg = "discrete",
#                       verbose = TRUE, null_diagnostic_plots = T)
#   end_discrete <- proc.time() - start_discrete
#   # 37, 30, 31, 33, 140, 33, 33, 51, 57, 90
# 
#   end_sand[3]; end_discrete[3]
#   # sandwich took 140 seconds, discrete took 534
# 
#   lik_sand <- sapply(corn_sand1_new$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
#   lik_discrete <- sapply(corn_discrete1$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
#   plot(lik_discrete - lik_sand, corn_sand1_new$coef$score_stat[1:10] - corn_discrete1$coef$score_stat[1:10])
# 
#   summary(lik_discrete - lik_sand)
#   summary(corn_sand1_new$coef$score_stat[1:30] - corn_discrete1$coef$score_stat[1:30])
# 
#   # n = 64, J = 479, p = 2
#   # sandwich takes 140 seconds, discrete takes 534 seconds
#   # ll is nearly always higher (up to 55) for discrete, score stats are similar, with
#   # one off by 0.8 and one off by 0.3
# })

# # commenting out to avoid note about corncob data not being part of package
# # however keeping because these are useful for testing
# test_that("test corncob data, n = 64, J = 242, p = 2", {
#   
#   skip("don't want to test automatically")
#   
#   # corncob data
#   
#   library(corncob)
#   library(phyloseq)
#   data(soil_phylo_sample)
#   data(soil_phylo_otu)
#   soil_phylo <- phyloseq::phyloseq(phyloseq::sample_data(soil_phylo_sample),
#                                    phyloseq::otu_table(soil_phylo_otu, taxa_are_rows = TRUE),
#                                    phyloseq::tax_table(soil_phylo_taxa))
#   soil_phylo <- subset_samples(soil_phylo, Day %in% 0:1 & Amdmt %in% 0:1)
#   soil_phylo_firm <- subset_taxa(soil_phylo, Phylum == "Firmicutes")
#   soil_samp <- data.frame(sample_data(soil_phylo_firm))
#   soil_otu <- t(unclass(otu_table(soil_phylo_firm)))
#   
#   to_rm <- which(colSums(soil_otu) == 0)
#   soil_otu <- soil_otu[, -to_rm]
#   
#   soil_fit <- emuFit(formula = ~ Amdmt,
#                      data = soil_samp,
#                      Y = soil_otu,
#                      run_score_tests = FALSE, compute_cis = FALSE,
#                      tolerance = 1e-5, verbose = "development")
#   start_sand <- proc.time()
#   corn_sand1_new <- emuFit(formula = ~ Amdmt,
#                            data = soil_samp,
#                            Y = soil_otu,
#                            fitted_model = soil_fit,
#                            refit = FALSE,
#                            compute_cis = FALSE,
#                            test_kj = data.frame(k = 2, j = 1:30), null_fit_alg = "constraint_sandwich",
#                            verbose = TRUE, null_diagnostic_plots = T)
#   end_sand <- proc.time() - start_sand
#   # 11, 2, 4, 4, 1, 4, 7, 3, 5, 2, 7, 2, 9, 2, 2, 3, 2, 3, 6, 2, 2, 3, 2, 2, 3, 3, 2, 2, 2, 3
#   start_discrete <- proc.time()
#   corn_discrete1 <- emuFit(formula = ~ Amdmt,
#                            data = soil_samp,
#                            Y = soil_otu,
#                            fitted_model = soil_fit,
#                            refit = FALSE,
#                            compute_cis = FALSE,
#                            test_kj = data.frame(k = 2, j = 1:30), null_fit_alg = "discrete",
#                            verbose = TRUE, null_diagnostic_plots = T)
#   end_discrete <- proc.time() - start_discrete
#   # 4, 2, 3, 2, 2, 3, 5, 2, 6, 2, 4, 2, 8, 2, 2, 2, 2, 3, 7, 2, 2, 2, 3, 3, 2, 2, 2, 3, 2, 4
#   
#   end_sand[3]; end_discrete[3]
#   # sandwich takes 101 seconds, discrete takes 95
#   
#   lik_sand <- sapply(corn_sand1_new$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
#   lik_discrete <- sapply(corn_discrete1$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
#   plot(lik_discrete - lik_sand, corn_sand1_new$coef$score_stat[1:30] - corn_discrete1$coef$score_stat[1:30])
#   
#   summary(lik_discrete - lik_sand)
#   summary(corn_sand1_new$coef$score_stat[1:30] - corn_discrete1$coef$score_stat[1:30])
#   
#   # n = 64, J = 242, p = 2
#   # sandwich takes 101 seconds, discrete takes 95 seconds
#   # ll is quite similar for both methods (-1.5 -> 5), and score stats are mostly similar with 3 that
#   # vary more (up to magnitude of 1.2 difference)
# })

# # commenting out to avoid note about corncob data not being part of package
# # however keeping because these are useful for testing
# test_that("test corncob data, n = 119, J = 264, p = 3", {
#   
#   skip("don't want to test automatically")
#   
#   # corncob data
#   
#   library(corncob)
#   library(phyloseq)
#   data(soil_phylo_sample)
#   data(soil_phylo_otu)
#   soil_phylo <- phyloseq::phyloseq(phyloseq::sample_data(soil_phylo_sample),
#                                    phyloseq::otu_table(soil_phylo_otu, taxa_are_rows = TRUE),
#                                    phyloseq::tax_table(soil_phylo_taxa))
#   #soil_phylo <- subset_samples(soil_phylo, Day %in% 0:1 & Amdmt %in% 0:1)
#   soil_phylo_firm <- subset_taxa(soil_phylo, Phylum == "Firmicutes")
#   soil_samp <- data.frame(sample_data(soil_phylo_firm))
#   soil_otu <- t(unclass(otu_table(soil_phylo_firm)))
#   soil_otu <- soil_otu[, ]
#   
#   soil_fit <- emuFit(formula = ~ Amdmt,
#                      data = soil_samp,
#                      Y = soil_otu,
#                      run_score_tests = FALSE, compute_cis = FALSE,
#                      tolerance = 1e-5, verbose = "development")
#   start_sand <- proc.time()
#   corn_sand1_new <- emuFit(formula = ~ Amdmt,
#                            data = soil_samp,
#                            Y = soil_otu,
#                            fitted_model = soil_fit,
#                            refit = FALSE,
#                            compute_cis = FALSE,
#                            test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "constraint_sandwich",
#                            verbose = TRUE, null_diagnostic_plots = T)
#   end_sand <- proc.time() - start_sand
#   # 4, 7, 6, 7, 8, 5, 10, 5, 8, 3, 9, 5, 10, 7, 4, 6, 5, 4, 9, 8, 
#   start_discrete <- proc.time()
#   corn_discrete1 <- emuFit(formula = ~ Amdmt,
#                            data = soil_samp,
#                            Y = soil_otu,
#                            fitted_model = soil_fit,
#                            refit = FALSE,
#                            compute_cis = FALSE,
#                            test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "discrete",
#                            verbose = TRUE, null_diagnostic_plots = T)
#   end_discrete <- proc.time() - start_discrete
#   # 13, 13, 13, 12, 14, 12, 23, 12, 23, 12, 25, 13, 25, 17, 12, 13, 13, 13, 29, 12
#   
#   end_sand[3]; end_discrete[3]
#   # sandwich takes 128 seconds, discrete takes 318
#   
#   lik_sand <- sapply(corn_sand1_new$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
#   lik_discrete <- sapply(corn_discrete1$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
#   plot(lik_discrete - lik_sand, corn_sand1_new$coef$score_stat[1:20] - corn_discrete1$coef$score_stat[1:20])
#   
#   summary(lik_discrete - lik_sand)
#   summary(corn_sand1_new$coef$score_stat[1:20] - corn_discrete1$coef$score_stat[1:20])
#   
#   # n = 119, J = 264, p = 3
#   # sandwich takes 128 seconds, discrete takes 318 seconds
#   # ll is quite similar for both methods (discrete up to 2 higher), and score stats are mostly similar
# })
# 
# # commenting out to avoid note about corncob data not being part of package
# # however keeping because these are useful for testing
# test_that("test corncob data, n = 119, J = 109, p = 3", {
#   
#   skip("don't want to test automatically")
#   
#   # corncob data
#   
#   library(corncob)
#   library(phyloseq)
#   data(soil_phylo_sample)
#   data(soil_phylo_otu)
#   soil_phylo <- phyloseq::phyloseq(phyloseq::sample_data(soil_phylo_sample),
#                                    phyloseq::otu_table(soil_phylo_otu, taxa_are_rows = TRUE),
#                                    phyloseq::tax_table(soil_phylo_taxa))
#   #soil_phylo <- subset_samples(soil_phylo, Day %in% 0:1 & Amdmt %in% 0:1)
#   soil_phylo_elus <- subset_taxa(soil_phylo, Phylum == "Elusimicrobia")
#   soil_samp <- data.frame(sample_data(soil_phylo_elus))
#   soil_otu <- t(unclass(otu_table(soil_phylo_elus)))
#   soil_otu <- soil_otu[, ]
#   
#   soil_fit <- emuFit(formula = ~ Amdmt,
#                      data = soil_samp,
#                      Y = soil_otu,
#                      run_score_tests = FALSE, compute_cis = FALSE,
#                      tolerance = 1e-5, verbose = "development")
#   start_sand <- proc.time()
#   corn_sand1_new <- emuFit(formula = ~ Amdmt,
#                            data = soil_samp,
#                            Y = soil_otu,
#                            fitted_model = soil_fit,
#                            refit = FALSE,
#                            compute_cis = FALSE,
#                            test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "constraint_sandwich",
#                            verbose = TRUE, null_diagnostic_plots = T)
#   end_sand <- proc.time() - start_sand
#   # 3, 2, 2, 2, 2, 2, 2, 3, 1, 1, 3, 2, 1, 2, 2, 1, 1, 1, 1, 2
#   start_discrete <- proc.time()
#   corn_discrete1 <- emuFit(formula = ~ Amdmt,
#                            data = soil_samp,
#                            Y = soil_otu,
#                            fitted_model = soil_fit,
#                            refit = FALSE,
#                            compute_cis = FALSE,
#                            test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "discrete",
#                            verbose = TRUE, null_diagnostic_plots = T)
#   end_discrete <- proc.time() - start_discrete
#   # 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1
#   
#   end_sand[3]; end_discrete[3]
#   # sandwich takes 35 seconds, discrete takes 12
#   
#   lik_sand <- sapply(corn_sand1_new$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
#   lik_discrete <- sapply(corn_discrete1$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
#   plot(lik_discrete - lik_sand, corn_sand1_new$coef$score_stat[1:20] - corn_discrete1$coef$score_stat[1:20])
#   
#   summary(lik_discrete - lik_sand)
#   summary(corn_sand1_new$coef$score_stat[1:20] - corn_discrete1$coef$score_stat[1:20])
#   
#   # n = 119, J = 109, p = 3
#   # sandwich takes 35 seconds, discrete takes 12 seconds
#   # ll and score stats are very similar
# })

# # commenting out to avoid note about corncob data not being part of package
# # however keeping because these are useful for testing
# test_that("test corncob data, n = 119, J = 147, p = 3", {
#   
#   skip("don't want to test automatically")
#   
#   # corncob data
#   
#   library(corncob)
#   library(phyloseq)
#   data(soil_phylo_sample)
#   data(soil_phylo_otu)
#   soil_phylo <- phyloseq::phyloseq(phyloseq::sample_data(soil_phylo_sample),
#                                    phyloseq::otu_table(soil_phylo_otu, taxa_are_rows = TRUE),
#                                    phyloseq::tax_table(soil_phylo_taxa))
#   #soil_phylo <- subset_samples(soil_phylo, Day %in% 0:1 & Amdmt %in% 0:1)
#   soil_phylo_gemm <- subset_taxa(soil_phylo, Phylum == "Gemmatimonadetes")
#   soil_samp <- data.frame(sample_data(soil_phylo_gemm))
#   soil_otu <- t(unclass(otu_table(soil_phylo_gemm)))
#   soil_otu <- soil_otu[, ]
#   
#   soil_fit <- emuFit(formula = ~ Amdmt,
#                      data = soil_samp,
#                      Y = soil_otu,
#                      run_score_tests = FALSE, compute_cis = FALSE,
#                      tolerance = 1e-5, verbose = "development")
#   start_sand <- proc.time()
#   corn_sand1_new <- emuFit(formula = ~ Amdmt,
#                            data = soil_samp,
#                            Y = soil_otu,
#                            fitted_model = soil_fit,
#                            refit = FALSE,
#                            compute_cis = FALSE,
#                            test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "constraint_sandwich",
#                            verbose = TRUE, null_diagnostic_plots = T)
#   end_sand <- proc.time() - start_sand
#   # 3, 2, 3, 2, 3, 2, 2, 3, 4, 2, 4, 4, 3, 1, 1, 2, 3, 2, 2, 2
#   start_discrete <- proc.time()
#   corn_discrete1 <- emuFit(formula = ~ Amdmt,
#                            data = soil_samp,
#                            Y = soil_otu,
#                            fitted_model = soil_fit,
#                            refit = FALSE,
#                            compute_cis = FALSE,
#                            test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "discrete",
#                            verbose = TRUE, null_diagnostic_plots = T)
#   end_discrete <- proc.time() - start_discrete
#   # 3, 2, 3, 2, 2, 3, 2, 3, 3, 2, 5, 3, 2, 2, 2, 3, 3, 2, 2, 2
#   
#   end_sand[3]; end_discrete[3]
#   # sandwich takes 50 seconds, discrete takes 49
#   
#   lik_sand <- sapply(corn_sand1_new$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
#   lik_discrete <- sapply(corn_discrete1$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
#   plot(lik_discrete - lik_sand, corn_sand1_new$coef$score_stat[1:20] - corn_discrete1$coef$score_stat[1:20])
#   
#   summary(lik_discrete - lik_sand)
#   summary(corn_sand1_new$coef$score_stat[1:20] - corn_discrete1$coef$score_stat[1:20])
#   
#   # n = 119, J = 147, p = 3
#   # sandwich takes 50 seconds, discrete takes 49 seconds
#   # ll and score stats are very similar
# })