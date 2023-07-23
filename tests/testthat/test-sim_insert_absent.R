# test_that("We are inserting rows for missing taxa correctly", {
#
#   J <- 240
#   b1_magnitude <- 0
#   set.seed(0)
#   n_per_group <- 25
#   nboot <- 2
#   cor_structs <-
#     generate_sim_setting(J,
#                          n_correlations = 10,
#                          max_no_correlated = 25)
#
#   delta <- rnorm(J,sd = 2)
#
#   b1_pre <- rep(0,J) #obv adjust this for alternative
#
#   Y <- generate_sim_data(n_per_group = n_per_group,
#                          b0 = c(rep(-9,J/10),seq(-5,-2,length.out = 0.8*J),rep(2,J/10)),
#                          b1 = b1_pre*b1_magnitude,
#                          delta = delta,
#                          z_mean = 4,
#                          nb_size = .125,
#                          cor_structs = cor_structs)
#
#   absent <- apply(Y,2,function(x) sum(x) ==0)
#   which_absent <- which(absent)
#   X <- cbind(1,rep(c(0,1),each = n_per_group))
#
#   Y_arch <- Y
#   message("Removing absent taxa")
#   Y <- Y[,!absent]
#
#   ml_fit <- emuFit(X = X, Y = Y,method = "ML")
#
#
#   ml_cis <- emuCI(emuMod = ml_fit,nboot = nboot,
#                   ninner = 10)
#   ml_cis <- ml_cis[ml_cis$row ==2,]
#
#   ml_fit_adj <- emuFit(X = X, Y = Y_arch + 0.01,
#                        method = "ML")
#
#   ml_adj_cis <- emuCI(ml_fit_adj,nboot = nboot)
#
#   ml_adj_cis <- ml_adj_cis[ml_adj_cis$row==2,]
#
#   ml_cis_complete <- sim_insert_absent(ci_object = ml_cis,
#                               which_absent = which_absent,
#                               J = J)
#
#   expect_equal(which_absent,which(is.na(ml_cis_complete$estimate)))
#
#   Y_in_both_groups <- apply(Y_arch[1:n_per_group,],2,function(x) sum(x>0)>0)&
#     apply(Y_arch[(n_per_group + 1):(2*n_per_group),],2,function(x) sum(x>0)>0)
#   expect_true(max(abs(
#     ml_cis_complete$estimate - ml_adj_cis$estimate)[Y_in_both_groups])
#    <0.5)
# })
