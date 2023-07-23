# test_that("Empirical means in simulated model are consistent with theoretical values", {
#
#   J <- 60
#   n_per_group <- 50
#   set.seed(0)
#   cor_structs <-
#     generate_sim_setting(J,
#                        n_correlations = 20,
#                        max_no_correlated = 10)
#
#   delta <- rnorm(J,sd = 3)
#
#
#  draws <- lapply(1:100,
#         function(k)
#           generate_sim_data(n_per_group = 50,
#                     b0 = c(rep(-9,J/10),seq(-5,-2,length.out = 0.8*J),rep(2,J/10)),
#                     b1 = rep(0,J),
#                     delta = delta,
#                     z_mean = 4,
#                     z_sd = 2,
#                     nb_size = 1,
#                     cor_structs = cor_structs,
#                     constant_z = TRUE))
#
#   means <- matrix(nrow = 2*n_per_group,
#                   ncol = J)
#  for(i in 1:(2*n_per_group)){
#    for(j in 1:J){
#      means[i,j] <- mean(sapply(1:100,
#                                function(k) draws[[k]][i,j]))
#    }
#  }
#
#   expect_true(abs(lm(log(apply(means,2,mean))[-(1:3)]~delta[-(1:3)])$coef[2] - 1)< 0.05)
# })
#
# #
# test_that("Under null simulation settings, observed counts are reasonable", {
#
#   J <- 60
#   n_per_group <- 25
#   set.seed(0)
#   cor_structs <-
#     generate_sim_setting(J,
#                          n_correlations = 20,
#                          max_no_correlated = 10)
#
#   delta <- rnorm(J,sd = 2)
#
#   sparsity <- numeric(10)
#   row_sparsity <- numeric(10)
#   maxes <- numeric(10)
#   sums <- numeric(10)
#   unobs <- numeric(10)
#   b1_pre <- rnorm(J)
#   b1_magnitude <- 0
#   for(draw in 1:10){
#   Y <- generate_sim_data(n_per_group = n_per_group,
#                          b0 = c(rep(-9,J/10),seq(-5,-2,length.out = 0.8*J),rep(2,J/10)),
#                          b1 = b1_pre*b1_magnitude,
#                          delta = delta,
#                          z_mean = 4,
#                          nb_size = .125,
#                          cor_structs = cor_structs)
#
#   Y_circ <- generate_sim_data(n_per_group = 25,
#                               b0 = c(rep(-9,J/10),seq(-5,-2,length.out = 0.8*J),rep(2,J/10)),
#                               b1 = rep(0,J),
#                               delta = delta,
#                               z_mean = 4,
#                               nb_size = .125,
#                               cor_structs = cor_structs,
#                               return_Y_circ= TRUE)
#
#   colnames(Y_circ) <- paste("taxon",1:60,sep = "")
#   rownames(Y_circ) <- paste("sample",1:50,sep = "")
#
#   # Y_circ %>%
#   #   as.data.frame() %>%
#   #   tibble::rownames_to_column("sample") %>%
#   #   mutate(group = rep(c("A","B"), each = 25)) %>%
#   #   tidyr::pivot_longer(-c(sample,group)) %>%
#   #   group_by(group,name) %>%
#   #   summarize(value = mean(value)) %>%
#   #   ungroup() %>%
#   #   pivot_wider(names_from = group, values_from = value) %>%
#   #   ggplot() +
#   #   geom_jitter(aes(x = A, y = B, color = name),
#   #               height = 0, width = .2) +
#   #   scale_y_log10() +
#   #   scale_x_log10() +
#   #   geom_abline(aes(slope = 1, intercept = 0),linetype = 2)
#
#   row_sparsity[draw] <- mean(apply(Y,2,sum)==0)
#   sparsity[draw] <- mean(Y ==0)
#   maxes[draw] <- max(Y)
#   sums[draw] <- sum(Y)
#   unobs[draw] <- sum(apply(Y, 2, function(x) mean(x ==0))==1)
#   }
#
#   expect_equal(mean(row_sparsity),0.07,tolerance = .1)
#   expect_equal(mean(sparsity),0.8,tolerance = 0.1)
#   expect_true(max(maxes)<1e7 & min(maxes)>.5e4)
#   expect_equal(mean(unobs),4.5,tolerance = 1)
#
#   J <- 250
#   n_per_group <- 25
#   set.seed(0)
#   cor_structs <-
#     generate_sim_setting(J,
#                          n_correlations = 20,
#                          max_no_correlated = 10)
#
#   delta <- rnorm(J,sd = 2)
#   for(draw in 1:10){
#     Y <- generate_sim_data(n_per_group = 25,
#                            b0 = c(rep(-9,J/10),seq(-5,-2,length.out = 0.8*J),rep(2,J/10)),
#                            b1 = rep(0,J),
#                            delta = delta,
#                            z_mean = 4,
#                            nb_size = .125,
#                            cor_structs = cor_structs,
#                            return_Y_circ= FALSE)
#
#     row_sparsity[draw] <- mean(apply(Y,2,sum)==0)
#     sparsity[draw] <- mean(Y ==0)
#     row_sparsity[draw] <- mean(apply(Y,2,sum)==0)
#     maxes[draw] <- max(Y)
#     sums[draw] <- sum(Y)
#     unobs[draw] <- sum(apply(Y, 2, function(x) mean(x ==0))==1)
#   }
#
#   expect_equal(mean(row_sparsity),0.073,tolerance= 0.1)
#   expect_equal(mean(sparsity),0.8,tolerance = 0.101)
#   expect_true(max(maxes)<1e7 & min(maxes)>1e4)
#   expect_equal(mean(unobs)/250,0.07,tolerance = 0.1)
#
#
# })
#
#
# test_that("We are actually inducing correlation", {
#
#   J <- 50
#   n_per_group <- 50
#   set.seed(0)
#   cor_structs <-
#     generate_sim_setting(J,
#                          n_correlations = 10,
#                          max_no_correlated = 10)
#
#   delta <- rnorm(J,sd = 2)
#
#   sparsity <- numeric(10)
#   maxes <- numeric(10)
#   sums <- numeric(10)
#   unobs <- numeric(10)
#     Y_circ <- generate_sim_data(n_per_group = 500,
#                                 b0 = c(rep(-9,J/10),seq(-5,-2,length.out = 0.8*J),rep(2,J/10)),
#                                 b1 = rep(0,J),
#                                 delta = delta,
#                                 z_mean = 4,
#                                 nb_size = .125,
#                                 cor_structs = cor_structs,
#                                 return_Y_circ= TRUE)
#
#     expect_equal(min(cor(Y_circ)),-.79,tolerance = 0.05)
#     expect_equal(max(cor(Y_circ)[upper.tri(cor(Y_circ))]),0.79,
#                  tolerance = 0.05)
#
# })
#
# # test_that("Null simulations for choices of n = 25 and J = 50 work", {
# #
# #   J <- 60
# #   set.seed(0)
# #   cor_structs <-
# #     generate_sim_setting(J,
# #                          n_correlations = 10,
# #                          max_no_correlated = 10)
# #
# #   delta <- rnorm(J,sd = 2)
# #
# #
# #   n25_J60 <-
# #     do_one_sim(n_per_group = 25,
# #                          J = 60,
# #                          delta = delta,
# #                          b1_magnitude = 0,
# #                          cor_structs = cor_structs,
# #                seed = 0,
# #                          nboot = 2)
# #
# #
# #   n100_J60 <-
# #     do_one_sim(n_per_group = 100,
# #                J = 60,
# #                delta = delta,
# #                b1_magnitude = 0,
# #                cor_structs = cor_structs,
# #                seed = 0,
# #                nboot = 2)
# #
# #   J <- 240
# #   set.seed(0)
# #   cor_structs <-
# #     generate_sim_setting(J,
# #                          n_correlations = 10,
# #                          max_no_correlated = 25)
# #
# #   delta <- rnorm(J,sd = 2)
# #
# #   n25_J240 <-
# #     do_one_sim(n_per_group = 25,
# #                J = 240,
# #                delta = delta,
# #                b1_magnitude = 0,
# #                cor_structs = cor_structs,
# #                seed = 0,
# #                nboot = 2)
# #
# #   n100_J240 <-
# #     do_one_sim(n_per_group = 100,
# #                J = 240,
# #                delta = delta,
# #                b1_magnitude = 0,
# #                cor_structs = cor_structs,
# #                seed = 0,
# #                nboot = 2)
# #   expect_true(nrow(n25_J60)== 60*4)
# #   expect_true(nrow(n100_J60) == 60*4)
# #   expect_true(nrow(n25_J240)== 240*4)
# #   expect_true(nrow(n100_J240) == 240*4)
# #
# # })
