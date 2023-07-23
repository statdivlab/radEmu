#
# test_that("emuCI with type = `sandwich` returns reasonable output with beta_J = 0 constraint",{
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20))
#   z <- rnorm(40) +8
#   b0 <- rnorm(10)
#   b1 <- 1:10 - 5.5
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 10, nrow = 40)
#
#   sim_results <- vector("list",100)
#   # for(sim in 1:100){
#   #   print(sim)
#   # set.seed(sim)
#   for(i in 1:40){
#     for(j in 1:10){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- #rpois(1, lambda = temp_mean)
#         rnbinom(1, mu = temp_mean, size=.5)
#     }
#   }
#
#   fit_fl <- emuFit_clone(Y = Y, X = X,
#                          # tolerance = 1e-4,
#                constraint_fn = function(x) x[10])
#
#   fit_fl_cis <- emuCI(emuMod = fit_fl,
#                       conf_level = 0.95,
#                       verbose = TRUE,
#                       type = "sandwich")
#
#   # fit_fl_cis_boot <- emuCI(emuMod = fit_fl,
#   #                          conf_level = 0.95,
#   #                          nboot = 100,
#   #                          type = "iterated")
#
#
#   expect_equal(mean(fit_fl_cis$lower[10:18]<= -(9:1)&
#                       fit_fl_cis$upper[10:18]>= -(9:1 )),
#                0.88,tolerance = .1)
#
# })
#
#
#
#
#
#
# test_that("emuCI returns appropriately formatted output & is basically reasonable
# on negative binomial test data", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 50))
#   z <- rnorm(100) +8
#   b0 <- rnorm(10)
#   b1 <- 1:10
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 10, nrow = 100)
#
#   sim_results <- vector("list",100)
#   # for(sim in 1:100){
#   #   print(sim)
#     # set.seed(sim)
#   for(i in 1:100){
#     for(j in 1:10){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- #rpois(1, lambda = temp_mean)
#         rnbinom(1, mu = temp_mean, size=0.25)
#     }
#   }
#   # constraint_fn <- (function(x){ x[1]})
#   fl_fit <- emuFit(X = X, Y = Y,method = "FL",
#                    fast_fit = TRUE,
#                    B = matrix(0,ncol = 10, nrow = 2),
#                    # constraint_fn = constraint_fn,
#                    reweight = FALSE,
#                                      tolerance = 0.01)
# #
# # fl_fit <- emuFit_clone(X = X, Y = Y,method = "FL",
# #                  # fast_fit = TRUE,
# #                  B = matrix(0,ncol = 10, nrow = 2),
# #                  constraint_fn = function(x) x[10],
# #                  # reweight = TRUE,
# #                  tolerance = 0.01)
#
#   # fl_fit <- emuFit_clone(X = X,
#   #                        Y = Y,
#   #                        method = "FL",
#   #                        constraint = function(x) x[10]
#   #                        )
#
#   boot_ci <-
#     emuCI(fl_fit,
#                    nboot = 25,
#                    ninner = 5,
#                    # ninner = 20,
#           seed = 2,
#                    type = "SBstudentized",
#           subsample_exp = sqrt(1/2)
#
#                    # parallel = TRUE,
#                    # ncores = 5
#   )
# #
# #   boot_ci %>%
# #     filter(row == 2) %>%
# #     ggplot() +
# #     geom_errorbar(aes(x = outcome_index, ymin = lower, ymax = upper),
# #                   width = 0.1) +
# #     geom_abline(aes(intercept = -5.5, slope = 1))
#
#
#   true_vals <- 1:10 - median(1:10)
#
#   # lapply(1:100, function(k) sim_results[[k]]$lower[11:20] <= true_vals &
#   #          sim_results[[k]]$upper[11:20]>= true_vals) %>%
#   #   do.call(rbind,.) %>% apply(2,mean)
#   expect_true(nrow(boot_ci)==20)
#   expect_true(ncol(boot_ci) == 7)
#   truevals <- b1 - median(b1)
#   expect_equal(mean((boot_ci$lower[11:20]<=truevals)&(boot_ci$upper[11:20]>=truevals)),0.85,
#                tolerance = .1)
# })
#
#
# test_that("emuCI returns appropriately formatted output & is basically reasonable
# on negative binomial test data (with reweight = FALSE)", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20))
#   z <- rnorm(40) +8
#   b0 <- rnorm(10)
#   b1 <- 1:10
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 10, nrow = 40)
#
#   for(i in 1:40){
#     for(j in 1:10){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- #rpois(1, lambda = temp_mean)
#         rnbinom(1, mu = temp_mean, size=0.25)
#     }
#   }
#   # constraint_fn <- (function(x){ x[1]})
#   fl_fit <- emuFit_clone(X = X, Y = Y,method = "FL",
#                    # fast_fit = TRUE,
#                    # constraint_fn = constraint_fn,
#                    reweight = FALSE,
#                    tolerance = 0.01)
#
#   boot_ci <- emuCI(fl_fit,
#                    nboot = 100,
#                    ninner = 5,
#                    # ninner = 20,
#                    type = "SBstudentized"
#                    # parallel = TRUE,
#                    # ncores = 5
#   )
#   # boot_ci %>%
#   #   filter(row ==1) %>%
#   #   ggplot() +
#   #   geom_errorbar(aes(x = outcome_index, ymin = lower, ymax = upper),
#   #                 width=0.25) +
#   #   geom_point(aes(x = outcome_index, y = estimate)) +
#   #   geom_abline(aes(slope=1,intercept=-5.5))
#   expect_true(nrow(boot_ci)==20)
#   expect_true(ncol(boot_ci) == 7)
#   truevals <-b1 - median(b1)
#   expect_equal(mean((boot_ci$lower[11:20]<=truevals)&(boot_ci$upper[11:20]>=truevals)),0.9,
#                tolerance = .2)
# })
#
#
# # test_that("emuCI works with fast_fit = TRUE", {
# #   set.seed(4323)
# #   X <- cbind(1,rep(c(0,1),each = 20))
# #   z <- rnorm(40) +8
# #   b0 <- rnorm(10)
# #   b1 <- 1:10
# #   b <- rbind(b0,b1)
# #   Y <- matrix(NA,ncol = 10, nrow = 40)
# #
# #   for(i in 1:40){
# #     for(j in 1:10){
# #       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
# #       Y[i,j] <- rpois(1, lambda = temp_mean)
# #     }
# #   }
# #   ml_fit <- fitted_model_ml<- emuFit(X = X, Y = Y,method = "ML",
# #                                      fast_fit = TRUE,
# #                                      tolerance = 0.01)
# #
# #   boot_ci <- emuCI(ml_fit,
# #                    nboot = 100,
# #                    fast_fit = TRUE,
# #                    rows_of_interest = c(1,2),
# #                    parallel = FALSE)
# #   expect_true(nrow(boot_ci)==20)
# #   expect_true(ncol(boot_ci) == 7)
# #   expect_true(boot_ci[20,"lower"]<4.5 & boot_ci[20,"upper"]>4.5)
# # })
#
# # test_that("emuCI works with fast_fit = TRUE and fix_augmentation = TRUE", {
# #   set.seed(4323)
# #   X <- cbind(1,rep(c(0,1),each = 20))
# #   z <- rnorm(40) +8
# #   b0 <- rnorm(10)
# #   b1 <- 1:10
# #   b <- rbind(b0,b1)
# #   Y <- matrix(NA,ncol = 10, nrow = 40)
# #
# #   for(i in 1:40){
# #     for(j in 1:10){
# #       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
# #       Y[i,j] <- rpois(1, lambda = temp_mean)
# #     }
# #   }
# #   fl_fit <- emuFit(X = X, Y = Y,
# #                    method = "FL",
# #                    fast_fit = TRUE,
# #                    tolerance = 0.01)
# #
# #   boot_ci <- emuCI(fl_fit,
# #                    nboot = 100,
# #                    fast_fit = TRUE,
# #                    fix_augmentation = TRUE,
# #                    rows_of_interest = c(1,2),
# #                    parallel = FALSE)
# #   expect_true(nrow(boot_ci)==20)
# #   expect_true(ncol(boot_ci) == 7)
# #   expect_true(boot_ci[20,"lower"]<4.5 & boot_ci[20,"upper"]>4.5)
# #   boot_ci$truth <- c(b0-median(b0),
# #                      b1-median(b1))
# #   expect_equal(mean((boot_ci$truth >= boot_ci$lower)&
# #     (boot_ci$truth <= boot_ci$upper)),0.95,tolerance = 0.1)
# # })
#
# # test_that("emuCI works with parallel = TRUE and ncore specified", {
# #   set.seed(4323)
# #   X <- cbind(1,rep(c(0,1),each = 20))
# #   z <- rnorm(40) +8
# #   b0 <- rnorm(10)
# #   b1 <- 1:10
# #   b <- rbind(b0,b1)
# #   Y <- matrix(NA,ncol = 10, nrow = 40)
# #
# #   for(i in 1:40){
# #     for(j in 1:10){
# #       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
# #       Y[i,j] <- rpois(1, lambda = temp_mean)
# #     }
# #   }
# #   ml_fit <- fitted_model_ml<- emuFit(X = X, Y = Y,method = "ML",
# #                                      tolerance = 0.01)
# #
# #   boot_ci <- emuCI(ml_fit,
# #                    nboot = 100,
# #                    rows_of_interest = c(1,2),
# #                    ncore = 5,
# #                    parallel = TRUE)
# #   expect_true(nrow(boot_ci)==20)
# #   expect_true(ncol(boot_ci) == 7)
# #   expect_true(boot_ci[20,"lower"]<4.5 & boot_ci[20,"upper"]>4.5)
# # })
