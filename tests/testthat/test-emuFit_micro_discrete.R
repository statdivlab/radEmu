
test_that("sim with single binary covariate", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  J <- 10
  n <- 40
  b1 <- 1:J - mean(1:J)
  Y <- radEmu:::simulate_data(n = n,
                              J = J,
                              X = X,
                              b0 = rnorm(J),
                              b1 = b1,
                              distn = "Poisson",
                              mean_z = 10)
  
  ml_fit <- emuFit_micro(X = X,
                         Y = Y,
                         constraint_fn = rep(list(function(x) mean(x)), 2), 
                         maxit = 200,
                         tolerance = 1e-3,
                         verbose = FALSE)
  ml_fit_fisher <- emuFit_micro(X = X,
                                Y = Y,
                                constraint_fn = rep(list(function(x) mean(x)), 2), 
                                maxit = 200,
                                tolerance = 1e-3,
                                verbose = FALSE,
                                use_discrete = FALSE)
  
  
  fit_discrete <- emuFit_micro_discrete(X = X,
                                        Y = Y,
                                        j_ref = 1)
  
  
  testthat::expect_lte(  max(abs(ml_fit -  
                                    rbind(fit_discrete[1,] - mean(fit_discrete[1,]), 
                                          fit_discrete[2,] - mean(fit_discrete[2,])))),
                          5e-3)
  
  testthat::expect_lte(max(abs(ml_fit - ml_fit_fisher)), 5e-3)
  
  ord <- sample(1:40, 40, replace = FALSE)
  X_reorder <- X[ord, ]
  Y_reorder <- Y[ord, ]
  
  ml_fit_reord <- emuFit_micro(X = X_reorder,
                               Y = Y_reorder,
                               constraint_fn = rep(list(function(x) mean(x)), 2), 
                               maxit = 200,
                               tolerance = 1e-3,
                               verbose = FALSE)
  ml_fit_fisher_reord <- emuFit_micro(X = X_reorder,
                                      Y = Y_reorder,
                                      constraint_fn = rep(list(function(x) mean(x)), 2), 
                                      maxit = 200,
                                      tolerance = 1e-3,
                                      verbose = FALSE,
                                      use_discrete = FALSE)
  
  
  fit_discrete_reord <- emuFit_micro_discrete(X = X_reorder,
                                              Y = Y_reorder,
                                              j_ref = 1)
  
  
  testthat::expect_lte(max(abs(fit_discrete - fit_discrete_reord)), 1e-8)
  testthat::expect_lte(max(abs(ml_fit - ml_fit_reord)), 1e-8)
  testthat::expect_lte(max(abs(ml_fit_reord - ml_fit_fisher_reord)), 5e-3)
  
})

test_that("sim with multiple-level covariate covariate", {
  set.seed(4323)
  dat <- data.frame(cat = rep(c("A", "B", "C"), each = 20))
  X <- make_design_matrix(formula = ~ cat, data = dat)
  X_reorder <- X[sample(1:60, 60), ]
  J <- 10
  n <- 60
  Y <- radEmu:::simulate_data(n = n,
                              J = J,
                              X = X,
                              B = rbind(rnorm(J, 0, 2), rnorm(J, 0, 2), rnorm(J, 0, 2)),
                              distn = "Poisson",
                              mean_z = 10)
  
  ml_fit <- emuFit_micro(X = X,
                         Y = Y,
                         constraint_fn = rep(list(function(x) mean(x)), 3), 
                         maxit = 200,
                         tolerance = 1e-3,
                         verbose = FALSE)
  ml_fit_fisher <- emuFit_micro(X = X,
                                Y = Y,
                                constraint_fn = rep(list(function(x) mean(x)), 3), 
                                maxit = 200,
                                tolerance = 1e-3,
                                verbose = FALSE,
                                use_discrete = FALSE)
  
  
  fit_discrete <- emuFit_micro_discrete(X = X,
                                        Y = Y,
                                        j_ref = 1)
  
  
  testthat::expect_lte(  max(abs(ml_fit -  
                                   rbind(fit_discrete[1,] - mean(fit_discrete[1,]), 
                                         fit_discrete[2,] - mean(fit_discrete[2,]),
                                         fit_discrete[3,] - mean(fit_discrete[3,])))),
                         5e-3)
  
  testthat::expect_lte(max(abs(ml_fit - ml_fit_fisher)), 5e-3)
  
  ord <- sample(1:60, 60, replace = FALSE)
  X_reorder <- X[ord, ]
  Y_reorder <- Y[ord, ]
  
  ml_fit_reord <- emuFit_micro(X = X_reorder,
                               Y = Y_reorder,
                               constraint_fn = rep(list(function(x) mean(x)), 3), 
                               maxit = 200,
                               tolerance = 1e-3,
                               verbose = FALSE)
  ml_fit_fisher_reord <- emuFit_micro(X = X_reorder,
                                      Y = Y_reorder,
                                      constraint_fn = rep(list(function(x) mean(x)), 3), 
                                      maxit = 200,
                                      tolerance = 1e-3,
                                      verbose = FALSE,
                                      use_discrete = FALSE)
  
  
  fit_discrete_reord <- emuFit_micro_discrete(X = X_reorder,
                                              Y = Y_reorder,
                                              j_ref = 1)
  
  
  testthat::expect_lte(max(abs(fit_discrete - fit_discrete_reord)), 1e-8)
  testthat::expect_lte(max(abs(ml_fit - ml_fit_reord)), 1e-8)
  testthat::expect_lte(max(abs(ml_fit_reord - ml_fit_fisher_reord)), 5e-3)
  
})

test_that("discrete works with wirbel data, binary design", {
  
  skip("too long for automated testing")
  
  data("wirbel_sample")
  data("wirbel_otu")
  data("wirbel_taxonomy")
  wirbel_phylo <- phyloseq::phyloseq(phyloseq::sample_data(wirbel_sample),
                                     phyloseq::otu_table(wirbel_otu, taxa_are_rows = FALSE),
                                     phyloseq::tax_table(wirbel_taxonomy))
  wirbel_genus <- phyloseq::tax_glom(wirbel_phylo, taxrank = "genus")
  zero_taxa <- phyloseq::taxa_sums(wirbel_genus) == 0
  wirbel_genus <- phyloseq::prune_taxa(zero_taxa == FALSE, wirbel_genus)
  wirbel_genus <- phyloseq::prune_samples(phyloseq::sample_sums(wirbel_genus) > 0, 
                                          wirbel_genus)
  design <- make_design_matrix(formula = ~ Group, Y = wirbel_genus)
  Y <- as.matrix(phyloseq::otu_table(wirbel_genus))
  
  ch_fit_old <- emuFit_micro_penalized(X = design, Y = Y,
                             constraint_fn = rep(list(pseudohuber_median), 6), 
                             use_discrete = FALSE)
  ch_fit <- emuFit_micro_penalized(X = design, Y = Y, 
                         constraint_fn = rep(list(pseudohuber_median), 6), 
                         use_discrete = TRUE)
  
  expect_true(all.equal(ch_fit_old$B, ch_fit$B, tol = 0.001))
  
})

test_that("discrete works with wirbel data, more categories", {
  
  skip("too long for automated testing")
  
  data("wirbel_sample")
  data("wirbel_otu")
  data("wirbel_taxonomy")
  wirbel_phylo <- phyloseq::phyloseq(phyloseq::sample_data(wirbel_sample),
                                     phyloseq::otu_table(wirbel_otu, taxa_are_rows = FALSE),
                                     phyloseq::tax_table(wirbel_taxonomy))
  wirbel_genus <- phyloseq::tax_glom(wirbel_phylo, taxrank = "genus")
  zero_taxa <- phyloseq::taxa_sums(wirbel_genus) == 0
  wirbel_genus <- phyloseq::prune_taxa(zero_taxa == FALSE, wirbel_genus)
  wirbel_genus <- phyloseq::prune_samples(phyloseq::sample_sums(wirbel_genus) > 0, 
                                          wirbel_genus)
  design <- make_design_matrix(formula = ~ Study, Y = wirbel_genus)
  Y <- as.matrix(phyloseq::otu_table(wirbel_genus))
  
  ch_fit_old <- emuFit_micro_penalized(X = design, Y = Y,
                                       constraint_fn = rep(list(pseudohuber_median), 6), 
                                       use_discrete = FALSE)
  ch_fit <- emuFit_micro_penalized(X = design, Y = Y, 
                                   constraint_fn = rep(list(pseudohuber_median), 6), 
                                   use_discrete = TRUE)
  
  expect_true(all.equal(ch_fit_old$B, ch_fit$B, tol = 0.001))
  
  ord <- sample(1:nrow(design), nrow(design), replace = FALSE)
  ch_fit_reord <- emuFit_micro_penalized(X = design[ord, ], Y = Y[ord, ], 
                                   constraint_fn = rep(list(pseudohuber_median), 6), 
                                   use_discrete = TRUE)
  
  expect_true(all.equal(ch_fit_reord$B, ch_fit$B, tol = 0.001))
  
})
