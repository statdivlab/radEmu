# test_that("multiplication works", {
#   J <- 60
#   n_per_group <- 25
#   nboot <- 2
#   set.seed(0)
#   verbose <- TRUE
#   cor_structs <-
#     generate_sim_setting(J,
#                                    n_correlations = 10,
#                                    max_no_correlated = 10
#
#   )
#
#   b1_pre <- c(1,rep(0,J/4 - 1),
#               -1,rep(0,J/4 - 1),
#               rep(0,J/4 - 1), 1,
#               rep(0,J/4 - 1),-1)
#
#   b1_magnitude <- 0
#
#   delta <- rnorm(J)
#
#
#
#
#
#   Y <- generate_sim_data(n_per_group = n_per_group,
#                          b0 = c(rep(-9,J/10),seq(-5,-2,length.out = 0.8*J),rep(2,J/10)),
#                          b1 = b1_pre*b1_magnitude,
#                          delta = delta,
#                          z_mean = 4,
#                          nb_size = 1,
#                          cor_structs = cor_structs)
#
#   absent <- apply(Y,2,function(x) sum(x) ==0)
#   which_absent <- which(absent)
#   X <- cbind(1,rep(c(0,1),each = n_per_group))
#
#   message("Removing absent taxa")
#   Y <- Y[,!absent]
#
#   constraint_fn <- (function(x){
#     median(x[(J/2 + 1):(J - sum(absent))])
#   })
#
#   ml_fit <- emuFit_clone(X = X, Y = Y,method = "ML",
#                    constraint_fn = constraint_fn,
#                    verbose = verbose)
#   ml_reweighted_fit <- emuFit_clone(X = X,
#                               Y = Y, method = "ML",
#                               B = ml_fit$B,
#                               constraint_fn = constraint_fn,
#                               reweight = TRUE,
#                               prefit = FALSE,
#                               verbose = verbose)
#
#   fl_fit <- emuFit_clone(X = X, Y = Y,
#                    method = "FL",
#                    # fast_fit = TRUE,
#                    constraint_fn = constraint_fn,
#                    verbose = verbose)
#
#   fl_reweighted_fit <- emuFit_clone(X = X,
#                               Y = Y,
#                               B = fl_fit$B,
#                               # fast_fit = TRUE,
#                               method = "FL",
#                               reweight = TRUE,
#                               prefit = FALSE,
#                               verbose = verbose)
#
#   ml_cis <- emuCI(ml_fit,nboot = nboot,
#                   # fast_fit = FALSE,
#                   ninner = 10,
#                   type = "SBstudentized",
#                   verbose = verbose)
#
#   ml_reweighted_cis <- emuCI(ml_reweighted_fit,nboot = nboot,
#                              type = "SBstudentized",
#                              ninner = 10,
#                              # fast_fit = FALSE,
#                              verbose = verbose)
#
#   fl_cis <- emuCI(fl_fit,nboot = nboot,
#                   type = "SBstudentized",
#                   ninner = 10,
#                   verbose = verbose)
#
#   fl_reweighted_cis <- emuCI(fl_reweighted_fit,nboot = nboot,
#                              ninner = 10,
#                              type = "SBstudentized",
#                              verbose = verbose)
#
#   n_missing <- length(which_absent)
#
#   message("Reinserting rows for absent taxa")
#
#   ml_cis <- sim_insert_absent(ci_object = ml_cis[ml_cis$row == 2,],
#                               which_absent = which_absent,
#                               J = J)
#
#   ml_reweighted_cis <- sim_insert_absent(
#     ci_object = ml_reweighted_cis[ml_reweighted_cis$row ==2,],
#                                          which_absent = which_absent,
#                                          J = J)
#
#   fl_cis <- sim_insert_absent(
#     ci_object = fl_cis[fl_cis$row ==2,],
#                               which_absent = which_absent,
#                               J = J)
#
#   fl_reweighted_cis <- sim_insert_absent(
#     ci_object = fl_reweighted_cis[fl_reweighted_cis$row ==2,],
#                                          which_absent = which_absent,
#                                          J = J)
#
#   message("Collating results")
#   ml_cis$estimator <- "ml"
#   ml_cis$weighted <- "no"
#   ml_reweighted_cis$estimator <- "ml"
#   ml_reweighted_cis$weighted <- "yes"
#   fl_cis$estimator <- "fl"
#   fl_cis$weighted <- "no"
#   fl_reweighted_cis$estimator <- "fl"
#   fl_reweighted_cis$weighted <- "yes"
#
#   # rbind(ml_cis,
#   #       ml_reweighted_cis,
#   #       fl_cis,
#   #       fl_reweighted_cis) %>%
#   #   ggplot() +
#   #   geom_errorbar(aes(x = outcome_index,
#   #                     ymin = lower,
#   #                     ymax= upper,
#   #                     group = interaction(weighted,estimator),
#   #                     color = estimator,
#   #                     linetype = weighted),
#   #                 position = position_dodge(1),
#   #                 width = .1) +
#   #   coord_cartesian(ylim = c(-5,5)) +
#   #   facet_grid(estimator~weighted) +
#   #   theme_bw()
#
#   rbind(ml_cis,
#         ml_reweighted_cis,
#         fl_cis,
#         fl_reweighted_cis) %>%
#     dplyr::filter(row == 2) %>%
#     dplyr::mutate(covers = lower<=0&upper>=0) %>%
#     dplyr::group_by(estimator,weighted) %>%
#     dplyr::summarize(coverage = mean(covers,na.rm = TRUE))
#
#   expect_true(nrow(  rbind(ml_cis,
#                            ml_reweighted_cis,
#                            fl_cis,
#                            fl_reweighted_cis) ) == 240)
#
#
# })
