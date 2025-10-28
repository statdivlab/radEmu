set.seed(59542234)
n <- 40
J <- 50
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
constraint_fn <- vector(mode = "list", length = 2)
constraint_grad_fn <- vector(mode = "list", length = 2)
reference_set <- list(1:30, 20:50)
for (k in 1:2) {
  constraint_fn[[k]] <- (function(k) {
    force(k)
    (function(x) { 
      pseudohuber_median(x[reference_set[[k]]])
    })})(k)
  constraint_grad_fn[[k]] <- (function(k) {
    force(k)
    (function(x) {
      grad <- rep(0, length(x))
      grad[reference_set[[k]]] <- dpseudohuber_median_dx(x[reference_set[[k]]])
      return(grad)
    })
  })(k)
}

test_that("we can choose null fit algorithm argument, subset", {
  fit1 <- emuFit(X = X, Y = Y, test_kj = data.frame(k = 2, j = 5), null_fit_alg = "constraint_sandwich", tolerance = 1e-2,
                 match_row_names = FALSE, constraint_fn = constraint_fn, constraint_grad_fn = constraint_grad_fn)
  fit2 <- emuFit(X = X, Y = Y, test_kj = data.frame(k = 2, j = 5), null_fit_alg = "augmented_lagrangian", tolerance = 1e-2,
                 match_row_names = FALSE, constraint_fn = constraint_fn, constraint_grad_fn = constraint_grad_fn)
  
  # expect pretty similar results
  expect_true(all.equal(fit1$coef$score_stat[5], fit2$coef$score_stat[5], tol = 0.05))
  # but not exactly the same, because different algorithms are being used
  expect_false(fit1$coef$score_stat[5] == fit2$coef$score_stat[5])
})

test_that("can run without errors", {
  
  skip("takes too long in automated testing")
  
  fit1 <- emuFit(X = X, Y = Y, test_kj = data.frame(k = 2, j = 1:50), null_fit_alg = "constraint_sandwich", tolerance = 1e-2,
                 match_row_names = FALSE, constraint_fn = constraint_fn, constraint_grad_fn = constraint_grad_fn,
                 verbose = TRUE, null_diagnostic_plots = TRUE)
  fit2 <- emuFit(X = X, Y = Y, test_kj = data.frame(k = 2, j = 1:50), null_fit_alg = "augmented_lagrangian", tolerance = 1e-2,
                 match_row_names = FALSE, constraint_fn = constraint_fn, constraint_grad_fn = constraint_grad_fn,
                 verbose = TRUE, null_diagnostic_plots = TRUE)
  
  fit1$coef$score_stat
  fit2$coef$score_stat
  
  plot((sapply(fit1$null_diagnostic_plots, function(x) tail(x$diagnostics_df$lik, 1)) - 
    sapply(fit2$null_diagnostic_plots, function(x) tail(x$diagnostics_df$lik, 1))), 
    (fit1$coef$score_stat - fit2$coef$score_stat))
  
  # result: constraint sandwich has higher ll than augmented lagrangian in 49/50 trials, 
  # max test stat difference is 0.12 and median differences is 0.009
  
})

# # commenting out to avoid note about corncob data not being part of package 
# # however keeping because these are useful for testing 
# 
# test_that("test corncob data", {
#   
#   skip("don't want to test automatically")
#   
#   # corncob data
# 
#   constraint_fn <- vector(mode = "list", length = 2)
#   constraint_grad_fn <- vector(mode = "list", length = 2)
#   reference_set <- list(1:20, 1:20, 30:50, 40:60)
#   for (k in 1:4) {
#     constraint_fn[[k]] <- (function(k) {
#       force(k)
#       (function(x) { 
#         pseudohuber_median(x[reference_set[[k]]])
#       })})(k)
#     constraint_grad_fn[[k]] <- (function(k) {
#       force(k)
#       (function(x) {
#         grad <- rep(0, length(x))
#         grad[reference_set[[k]]] <- dpseudohuber_median_dx(x[reference_set[[k]]])
#         return(grad)
#       })
#     })(k)
#   }
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
#                      tolerance = 1e-5, verbose = "development",
#                      constraint_fn = constraint_fn,
#                      constraint_grad_fn = constraint_grad_fn)
#   
#   corn_sand1_new <- emuFit(formula = ~ Amdmt + Day + Plants,
#                        data = soil_samp,
#                        Y = soil_otu,
#                        fitted_model = soil_fit,
#                        refit = FALSE,
#                        compute_cis = FALSE,
#                        test_kj = data.frame(k = 2, j = 1:30), null_fit_alg = "constraint_sandwich",
#                        verbose = TRUE, null_diagnostic_plots = T,
#                        constraint_fn = constraint_fn,
#                        constraint_grad_fn = constraint_grad_fn, maxit_null = 100)
#   
#   corn_aug1 <- emuFit(formula = ~ Amdmt + Day + Plants,
#                        data = soil_samp,
#                        Y = soil_otu,
#                        fitted_model = soil_fit,
#                        refit = FALSE,
#                        compute_cis = FALSE,
#                        test_kj = data.frame(k = 2, j = 1:30), null_fit_alg = "augmented_lagrangian",
#                        verbose = TRUE, null_diagnostic_plots = T,
#                       constraint_fn = constraint_fn,
#                       constraint_grad_fn = constraint_grad_fn)
# 
#   lik_sand <- sapply(corn_sand1_new$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
#   lik_aug <- sapply(corn_aug1$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
#   plot(lik_sand - lik_aug, corn_sand1_new$coef$score_stat[1:30] - corn_aug1$coef$score_stat[1:30])
# 
#   plot(corn_sand1_new$coef$score_stat[1:30], corn_sand1_new$coef$score_stat[1:30] - corn_aug1$coef$score_stat[1:30])
# 
#   sand_failed <- sapply(corn_sand1_new$null_diagnostic_plots, function(x) {"constraint_diff" %in% names(x$diagnostics_df)})
#   
#   summary(lik_aug - lik_sand)
#   summary((lik_aug - lik_sand)/lik_aug)
#   cbind(lik_aug - lik_sand, abs(corn_sand1_new$coef$score_stat[1:30] - corn_aug1$coef$score_stat[1:30]))
#   summary(abs(corn_sand1_new$coef$score_stat[1:30] - corn_aug1$coef$score_stat[1:30]))
#   
#   # result: constrained sandwich ll is lower than augmented lagrangian on average (although largest difference
#   # isn't very large), but test stats are different max 3.5, median 0.3
#   
#   # theory - the constraint sandwich alg doesn't work as well with such a small subset (size of 20),
#   # below tried with a larger subset size and the two algorithms give much more similar results
#   
# })

# test_that("test corncob data", {
# 
#   skip("don't want to test automatically")
# 
#   # corncob data
# 
#   constraint_fn <- vector(mode = "list", length = 2)
#   constraint_grad_fn <- vector(mode = "list", length = 2)
#   reference_set <- list(100:130, 1:30, 20:50, 40:70)
#   for (k in 1:4) {
#     constraint_fn[[k]] <- (function(k) {
#       force(k)
#       (function(x) {
#         pseudohuber_median(x[reference_set[[k]]])
#       })})(k)
#     constraint_grad_fn[[k]] <- (function(k) {
#       force(k)
#       (function(x) {
#         grad <- rep(0, length(x))
#         grad[reference_set[[k]]] <- dpseudohuber_median_dx(x[reference_set[[k]]])
#         return(grad)
#       })
#     })(k)
#   }
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
#                      tolerance = 1e-5, verbose = "development",
#                      constraint_fn = constraint_fn,
#                      constraint_grad_fn = constraint_grad_fn)
# 
#   corn_sand30 <- emuFit(formula = ~ Amdmt + Day + Plants,
#                        data = soil_samp,
#                        Y = soil_otu,
#                        fitted_model = soil_fit,
#                        refit = FALSE,
#                        compute_cis = FALSE,
#                        test_kj = data.frame(k = 2, j = 1:40), null_fit_alg = "constraint_sandwich",
#                        verbose = TRUE, null_diagnostic_plots = T,
#                        constraint_fn = constraint_fn,
#                        constraint_grad_fn = constraint_grad_fn, maxit_null = 100)
# 
#   corn_aug30 <- emuFit(formula = ~ Amdmt + Day + Plants,
#                        data = soil_samp,
#                        Y = soil_otu,
#                        fitted_model = soil_fit,
#                        refit = FALSE,
#                        compute_cis = FALSE,
#                        test_kj = data.frame(k = 2, j = 1:40), null_fit_alg = "augmented_lagrangian",
#                        verbose = TRUE, null_diagnostic_plots = T,
#                       constraint_fn = constraint_fn,
#                       constraint_grad_fn = constraint_grad_fn)
# 
#   lik_sand <- sapply(corn_sand30$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
#   lik_aug <- sapply(corn_aug30$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
#   plot(lik_sand - lik_aug, corn_sand30$coef$score_stat[1:40] - corn_aug30$coef$score_stat[1:40])
# 
#   summary(lik_aug - lik_sand)
#   summary((lik_aug - lik_sand)/lik_aug)
#   summary(abs(corn_sand30$coef$score_stat[1:40] - corn_aug30$coef$score_stat[1:40]))
# 
#   # result: constrained sandwich is ~4-5x faster
#   # ll differences are max 1
#   # test stats differ by median of 0.004, max of 0.1
# 
# })

# 
# test_that("test corncob data", {
#   
#   # result: constraint sandwich is faster, very similar likelihoods, similar test statistics 
#   
#   skip("don't want to test automatically")
#   
#   # corncob data
#   
#   constraint_fn <- vector(mode = "list", length = 2)
#   constraint_grad_fn <- vector(mode = "list", length = 2)
#   reference_set <- list(1:40, 1:40, 30:70, 40:80)
#   for (k in 1:4) {
#     constraint_fn[[k]] <- (function(k) {
#       force(k)
#       (function(x) { 
#         pseudohuber_median(x[reference_set[[k]]])
#       })})(k)
#     constraint_grad_fn[[k]] <- (function(k) {
#       force(k)
#       (function(x) {
#         grad <- rep(0, length(x))
#         grad[reference_set[[k]]] <- dpseudohuber_median_dx(x[reference_set[[k]]])
#         return(grad)
#       })
#     })(k)
#   }
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
#   soil_fit2 <- emuFit(formula = ~ Amdmt + Day + Plants,
#                      data = soil_samp,
#                      Y = soil_otu,
#                      run_score_tests = FALSE, compute_cis = FALSE,
#                      tolerance = 1e-5, verbose = "development",
#                      constraint_fn = constraint_fn,
#                      constraint_grad_fn = constraint_grad_fn)
#   
#   corn_sand2_new <- emuFit(formula = ~ Amdmt + Day + Plants,
#                            data = soil_samp,
#                            Y = soil_otu,
#                            fitted_model = soil_fit2,
#                            refit = FALSE,
#                            compute_cis = FALSE,
#                            test_kj = data.frame(k = 2, j = 1:30), null_fit_alg = "constraint_sandwich",
#                            verbose = TRUE, null_diagnostic_plots = T,
#                            constraint_fn = constraint_fn,
#                            constraint_grad_fn = constraint_grad_fn)
#   
#   corn_aug2 <- emuFit(formula = ~ Amdmt + Day + Plants,
#                       data = soil_samp,
#                       Y = soil_otu,
#                       fitted_model = soil_fit2,
#                       refit = FALSE,
#                       compute_cis = FALSE,
#                       test_kj = data.frame(k = 2, j = 1:30), null_fit_alg = "augmented_lagrangian",
#                       verbose = TRUE, null_diagnostic_plots = T,
#                       constraint_fn = constraint_fn,
#                       constraint_grad_fn = constraint_grad_fn)
#   
#   lik_sand <- sapply(corn_sand2_new$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
#   lik_aug <- sapply(corn_aug2$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
#   plot(lik_sand - lik_aug, abs(corn_sand2_new$coef$score_stat[1:30] - corn_aug2$coef$score_stat[1:30]))
#   plot(corn_sand2_new$coef$score_stat[1:30], abs(corn_sand2_new$coef$score_stat[1:30] - corn_aug2$coef$score_stat[1:30]))
#   
#   summary(lik_sand - lik_aug)
#   summary(abs(corn_sand2_new$coef$score_stat[1:30] - corn_aug2$coef$score_stat[1:30]))
#   
#   # results: ll from constrained sandwich - ll from augmented lagrangian ranges -0.5 - 0.06, median of 0.01, so quite similar
#   # test stats are max different by 0.06, median different by 0.003
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
  
  constraint_fn <- vector(mode = "list", length = 2)
  constraint_grad_fn <- vector(mode = "list", length = 2)
  reference_set <- lapply(1:6, function(x) {sample(1:ntaxa(wirbel_genus_ch), 40)})
  for (k in 1:6) {
    constraint_fn[[k]] <- (function(k) {
      force(k)
      (function(x) { 
        pseudohuber_median(x[reference_set[[k]]])
      })})(k)
    constraint_grad_fn[[k]] <- (function(k) {
      force(k)
      (function(x) {
        grad <- rep(0, length(x))
        grad[reference_set[[k]]] <- dpseudohuber_median_dx(x[reference_set[[k]]])
        return(grad)
      })
    })(k)
  }
  
  wirbel_full2 <- emuFit(formula = ~ Group + Gender + 
                          Age_spline.1 + Age_spline.2 + 
                          Sampling, 
                        Y = wirbel_genus_ch, run_score_tests = FALSE, 
                        compute_cis = FALSE, verbose = "development",
                        tolerance = 1e-5,
                        constraint_fn = constraint_fn,
                        constraint_grad_fn = constraint_grad_fn)
  
  start_sand <- proc.time()
  wirb_sand2 <- emuFit(formula = ~ Group + Gender + 
                        Age_spline.1 + Age_spline.2 + 
                        Sampling, 
                      Y = wirbel_genus_ch,
                      fitted_model = wirbel_full2,
                      refit = FALSE, 
                      compute_cis = FALSE, 
                      test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "constraint_sandwich",
                      verbose = FALSE, null_diagnostic_plots = T,
                      constraint_fn = constraint_fn,
                      constraint_grad_fn = constraint_grad_fn,
                      maxit_null = 100)
  end_sand <- proc.time() - start_sand
  # 
  start_aug <- proc.time()
  wirb_aug2 <- emuFit(formula = ~ Group + Gender + 
                       Age_spline.1 + Age_spline.2 + 
                       Sampling, 
                     Y = wirbel_genus_ch,
                     fitted_model = wirbel_full2,
                     refit = FALSE, 
                     compute_cis = FALSE, 
                     test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "augmented_lagrangian",
                     verbose = FALSE, null_diagnostic_plots = T,
                     constraint_fn = constraint_fn,
                     constraint_grad_fn = constraint_grad_fn)
  end_aug <- proc.time() - start_aug
  end_aug[3]; end_sand[3] 
  
  lik_sand <- sapply(wirb_sand2$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
  lik_aug <- sapply(wirb_aug2$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
  plot(lik_sand - lik_aug, wirb_sand2$coef$score_stat[1:20] - wirb_aug2$coef$score_stat[1:20])
  summary(lik_sand - lik_aug)
  summary(abs(wirb_sand2$coef$score_stat[1:20] - wirb_aug2$coef$score_stat[1:20]))
  
  # results: ll constrained sandwich - ll augmented lagrangian ranges -0.03 to 1.29, very close
  # test stats difference is max 0.25, median 0.01
  
  # uses random reference set of size 40
  
})

test_that("test wirbel data, reference set of size 30", {
  
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
  
  constraint_fn <- vector(mode = "list", length = 2)
  constraint_grad_fn <- vector(mode = "list", length = 2)
  reference_set <- lapply(1:6, function(x) {sample(1:ntaxa(wirbel_genus_ch), 30)})
  for (k in 1:6) {
    constraint_fn[[k]] <- (function(k) {
      force(k)
      (function(x) { 
        pseudohuber_median(x[reference_set[[k]]])
      })})(k)
    constraint_grad_fn[[k]] <- (function(k) {
      force(k)
      (function(x) {
        grad <- rep(0, length(x))
        grad[reference_set[[k]]] <- dpseudohuber_median_dx(x[reference_set[[k]]])
        return(grad)
      })
    })(k)
  }
  
  wirbel_full3 <- emuFit(formula = ~ Group + Gender + 
                           Age_spline.1 + Age_spline.2 + 
                           Sampling, 
                         Y = wirbel_genus_ch, run_score_tests = FALSE, 
                         compute_cis = FALSE, verbose = "development",
                         tolerance = 1e-5,
                         constraint_fn = constraint_fn,
                         constraint_grad_fn = constraint_grad_fn)
  
  start_sand <- proc.time()
  wirb_sand3 <- emuFit(formula = ~ Group + Gender + 
                         Age_spline.1 + Age_spline.2 + 
                         Sampling, 
                       Y = wirbel_genus_ch,
                       fitted_model = wirbel_full3,
                       refit = FALSE, 
                       compute_cis = FALSE, 
                       test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "constraint_sandwich",
                       verbose = TRUE, null_diagnostic_plots = T,
                       constraint_fn = constraint_fn,
                       constraint_grad_fn = constraint_grad_fn,
                       maxit_null = 100)
  end_sand <- proc.time() - start_sand
  # 
  start_aug <- proc.time()
  wirb_aug3 <- emuFit(formula = ~ Group + Gender + 
                        Age_spline.1 + Age_spline.2 + 
                        Sampling, 
                      Y = wirbel_genus_ch,
                      fitted_model = wirbel_full3,
                      refit = FALSE, 
                      compute_cis = FALSE, 
                      test_kj = data.frame(k = 2, j = 1:20), null_fit_alg = "augmented_lagrangian",
                      verbose = TRUE, null_diagnostic_plots = T,
                      constraint_fn = constraint_fn,
                      constraint_grad_fn = constraint_grad_fn)
  end_aug <- proc.time() - start_aug
  end_aug[3]; end_sand[3] 
  
  lik_sand <- sapply(wirb_sand3$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
  lik_aug <- sapply(wirb_aug3$null_diagnostic_plots, function(x) {tail(x$diagnostics_df$lik, 1)})
  plot(lik_sand - lik_aug, wirb_sand3$coef$score_stat[1:20] - wirb_aug3$coef$score_stat[1:20])
  summary(lik_sand - lik_aug)
  summary(abs(wirb_sand3$coef$score_stat[1:20] - wirb_aug3$coef$score_stat[1:20]))
  
  # results: constraint sandwich is ~4x faster
  # ll constrained sandwich - ll augmented lagrangian ranges -0.01 to 31
  # test stats difference is max 0.6, median 0.015
  # for larger differences, constrained sandwich ll is higher
  
  # uses random reference set of size 30
  
})

