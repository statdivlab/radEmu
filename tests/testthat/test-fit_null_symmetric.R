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
                 match_row_names = FALSE)
  fit2 <- emuFit(X = X, Y = Y, test_kj = data.frame(k = 2, j = 5), null_fit_alg = "augmented_lagrangian", tolerance = 1e-2,
                 match_row_names = FALSE)
  
  # expect pretty similar results
  expect_true(all.equal(fit1$coef$score_stat[5], fit2$coef$score_stat[5], tol = 0.01))
  # but not exactly the same, because different algorithms are being used
  expect_false(fit1$coef$score_stat[5] == fit2$coef$score_stat[5])
})

test_that("compare timing", {
  
  skip("long for automatic testing")
  
  sand_start <- proc.time()
  fit1 <- emuFit(X = X, Y = Y, test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "constraint_sandwich", tolerance = 1e-2,
                 match_row_names = FALSE, null_diagnostic_plots = T)
  sand_end <- proc.time() - sand_start
  lag_start <- proc.time()
  fit2 <- emuFit(X = X, Y = Y, test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "augmented_lagrangian", tolerance = 1e-2,
                 match_row_names = FALSE, null_diagnostic_plots = T)
  lag_end <- proc.time() - lag_start
  
  expect_true(lag_end[3] > sand_end[3])
  
  lik1 <- sapply(fit1$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
  lik2 <- sapply(fit2$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
  lik1 - lik2
  plot(lik1 - lik2, fit1$coef$score_stat - fit2$coef$score_stat)
  
  # when score stats are notably different, likelihood is higher for constraint_sandwich approach 
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
                 match_row_names = FALSE, null_diagnostic_plots = T)
  sand_end <- proc.time() - sand_start
  lag_start <- proc.time()
  fit2 <- emuFit(X = X, Y = Y_zinb, test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "augmented_lagrangian", tolerance = 1e-2,
                 match_row_names = FALSE, null_diagnostic_plots = T)
  lag_end <- proc.time() - lag_start
  
  expect_true(lag_end[3] > sand_end[3])
  lag_end[3] / sand_end[3] # augmented lagrangian is ~5x slower 
  
  lik1 <- sapply(fit1$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
  lik2 <- sapply(fit2$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
  summary(lik1 - lik2)
  summary((lik1 - lik2) / lik2)
  plot(lik1 - lik2, abs(fit1$coef$score_stat - fit2$coef$score_stat))
  summary(abs(fit1$coef$score_stat - fit2$coef$score_stat))
  # score stats are quite similar (median diff is 0.007, max diff is 0.08)
  # ll are quite similar, except for a few settings where constraint sandwich is lower by ~1000
  # this is still 1e-5 of overall ll
})

test_that("we get same null fit with different j_ref", {

  skip("confirmed this works but it is too slow for automated testing runs")
  
  k_constr <- 2
  j_constr <- 1
  p <- 2

  # constraint_fn <- rep(list(function(x){mean(x)}), 2)
  constraint_fn <- rep(list(function(x) pseudohuber_median(x, 0.1)), 2)
  # constraint_grad_fn <- function(x){dpseudohuber_median_dx(x,0.1)
  constraint_grad_fn <- rep(
    list(function(x) {
      dpseudohuber_median_dx(x, 0.1)
    }),
    2
  )

  full_fit <- emuFit_micro_penalized(
    X = X,
    Y = Y,
    B = NULL,
    constraint_fn = constraint_fn,
    tolerance = 1e-3,
    verbose = FALSE
  )

  B <- full_fit$B
  Y_aug <- full_fit$Y_augmented

  X_cup <- X_cup_from_X(X, J)

  j_ref <- 5
  null_fit <- fit_null(
    B = B,
    Y = Y_aug,
    X = X,
    X_cup = X_cup,
    k_constr = k_constr,
    j_constr = j_constr,
    j_ref = j_ref,
    constraint_fn = constraint_fn,
    constraint_grad_fn = constraint_grad_fn,
    constraint_tol = 1e-5,
    B_tol = 1e-4,
    verbose = FALSE,
    trackB = FALSE
  ) ## just track for one j

  null_repar_fit <- fit_null_symmetric(
    B = B,
    Y = Y_aug,
    X = X,
    X_cup = X_cup,
    j_constr = j_constr,
    k_constr = k_constr,
    j_ref = j_ref,
    constraint_fn = constraint_fn,
    constraint_grad_fn = constraint_grad_fn,
    B_tol = 1e-4,
    verbose = FALSE,
    use_optim = TRUE
  )

  null_repar_fit_fs <- fit_null_symmetric(
    B = B,
    Y = Y_aug,
    X = X,
    j_constr = j_constr,
    k_constr = k_constr,
    j_ref = j_ref,
    constraint_fn = constraint_fn,
    constraint_grad_fn = constraint_grad_fn,
    verbose = FALSE,
    use_optim = FALSE,
    X_cup = X_cup
  )

  #min_mse_lag finds the lagrange multiplier that minimizes the squared norm
  #of the derivative of the lagrangian (ll + lambda*(g - B[k_constr,j_constr])
  #and returns the minimized squared norm -- lower means more accurate fit under null

  null_min_lag_norm <- min_mse_lag(
    X = X,
    Y = Y,
    B = null_fit$B,
    constraint_grad_fn = constraint_grad_fn[[1]],
    k_constr = k_constr,
    j_constr = j_constr,
    j_ref = j_ref
  )
  null_repar_min_lag_norm <- min_mse_lag(
    X = X,
    Y = Y,
    B = null_repar_fit$B,
    constraint_grad_fn = constraint_grad_fn[[1]],
    k_constr = k_constr,
    j_constr = j_constr,
    j_ref = j_ref
  )
  null_repar_min_lag_norm_fs <- min_mse_lag(
    X = X,
    Y = Y,
    B = null_repar_fit_fs$B,
    constraint_grad_fn = constraint_grad_fn[[1]],
    k_constr = k_constr,
    j_constr = j_constr,
    j_ref = j_ref
  )

  #solns are at least close to equal
  # expect_equal(null_repar_fit$B, null_fit$B, tolerance = 1e-2)
  expect_equal(null_repar_fit_fs$B, null_fit$B, tolerance = 1e-2)

  #and to extent that they are not equal, repar fit is more accurate
  # expect_true(null_min_lag_norm > null_repar_min_lag_norm)
  expect_true(null_min_lag_norm > null_repar_min_lag_norm_fs)
})

# # commenting out to avoid note about corncob data not being part of package 
# # however keeping because these are useful for testing 
# test_that("test corncob data", {
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
#   soil_fit <- emuFit(formula = ~ Amdmt + Day + Plants,
#                      data = soil_samp,
#                      Y = soil_otu,
#                      run_score_tests = FALSE, compute_cis = FALSE,
#                      tolerance = 1e-5, verbose = "development")
#   start_sand <- proc.time()
#   corn_sand1_new <- emuFit(formula = ~ Amdmt + Day + Plants,
#                        data = soil_samp,
#                        Y = soil_otu,
#                        fitted_model = soil_fit,
#                        refit = FALSE,
#                        compute_cis = FALSE,
#                        test_kj = data.frame(k = 2, j = 1:30), null_fit_alg = "constraint_sandwich",
#                        verbose = TRUE, null_diagnostic_plots = T)
#   end_sand <- proc.time() - start_sand
#   # 8, 8, 6, 4, 7, 6, 7, 8, 9, 5, 5, 13, 9, 5, 8, 4, 3, 10, 4, 5, 6, 9, 8, 6, 5, 17, 4, 7, 6, 6
#   start_aug <- proc.time()
#   corn_aug1 <- emuFit(formula = ~ Amdmt + Day + Plants,
#                        data = soil_samp,
#                        Y = soil_otu,
#                        fitted_model = soil_fit,
#                        refit = FALSE,
#                        compute_cis = FALSE,
#                        test_kj = data.frame(k = 2, j = 1:30), null_fit_alg = "augmented_lagrangian",
#                        verbose = FALSE, null_diagnostic_plots = T)
#   end_aug <- proc.time() - start_aug
#   # 11, 91, 15, 11, 15, 16, 17, 90, 64, 16, 54, 140, 103, 15, 18, 13, 11, 18, 18, 12, 14, 17, 18, 20, 11, 33, 10, 16, 17, 18
#   end_aug[3]; end_sand[3]
# 
#   lik_sand <- sapply(corn_sand1_new$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
#   lik_aug <- sapply(corn_aug1$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
#   plot(lik_sand - lik_aug, corn_sand1_new$coef$score_stat[1:30] - corn_aug1$coef$score_stat[1:30])
# 
#   plot(corn_sand1_new$coef$score_stat[1:30], corn_sand1_new$coef$score_stat[1:30] - corn_aug1$coef$score_stat[1:30])
# 
#   summary(lik_aug - lik_sand)
#   summary(abs(corn_aug1$coef$score_stat - corn_sand1_new$coef$score_stat))
#   
#   # sandwich approach is ~3.5 times faster
#   # ll from two approaches ranges up to 3.5, overall very similar (for large scale of ll)
#   # test stats diff is median 0.015, max 0.28
#   # for the larger test stat differences, sandwich ll is larger 
#   
# })

test_that("test wirbel data", {
  
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
  ch_fit_full <- emuFit(formula = ~ Group + Gender + 
                          Age_spline.1 + Age_spline.2 + 
                          Sampling, 
                        Y = wirbel_genus_ch, run_score_tests = FALSE, 
                        compute_cis = FALSE, verbose = "development",
                        tolerance = 1e-5)
  
  start_sand <- proc.time()
  wirb_sand <- emuFit(formula = ~ Group + Gender + 
                        Age_spline.1 + Age_spline.2 + 
                        Sampling, 
                      Y = wirbel_genus_ch,
                      fitted_model = ch_fit_full,
                      refit = FALSE, 
                      compute_cis = FALSE, 
                      test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "constraint_sandwich",
                      verbose = TRUE, null_diagnostic_plots = T)
  end_sand <- proc.time() - start_sand
  # 8, 4, 4, 3, 6, 2, 7, 2, 4, 6, 3, 9, 6, 3, 4, 8, 2, 5, 5, 8
  start_aug <- proc.time()
  wirb_aug <- emuFit(formula = ~ Group + Gender + 
                       Age_spline.1 + Age_spline.2 + 
                       Sampling, 
                     Y = wirbel_genus_ch,
                     fitted_model = ch_fit_full,
                     refit = FALSE, 
                     compute_cis = FALSE, 
                     test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "augmented_lagrangian",
                     verbose = TRUE, null_diagnostic_plots = T)
  end_aug <- proc.time() - start_aug
  # 14, 11, 16, 12, 81, 12, 51, 14, 14, 93, 14, 10, 13, 14, 11, 47, 15, 35, 45, 55
  end_aug[3]; end_sand[3] # constraint sandwich is ~5.5x faster than augmented lagrangian
  
  lik_sand <- sapply(wirb_sand$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
  lik_aug <- sapply(wirb_aug$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
  plot(lik_sand - lik_aug, wirb_sand$coef$score_stat[1:20] - wirb_aug$coef$score_stat[1:20])
  
  summary(lik_sand - lik_aug)
  summary(abs(wirb_sand$coef$score_stat[1:20] - wirb_aug$coef$score_stat[1:20]))

  # results: constraint sandwich is ~5.5 faster
  # ll range up to ~45, constraint sandwich ll is typically larger (especially for larger differences)
  # test stat median difference is 0.01, max difference is 0.05
    
})
