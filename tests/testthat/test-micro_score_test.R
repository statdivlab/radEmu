
# #
# test_that("Test score test with Poisson data under null, J > n", {
#   set.seed(4323)
#   n <- 100
#   J <- 12
#   X <- cbind(1,rep(c(0,1),each = n/2))
#   z <- rnorm(n) + 8
#   b0 <- rnorm(J)
#   b1 <- seq(-1.5,1.5,length.out = J)
#   b1[J/2 + 0:1] <- 0
#   b <- rbind(b0,b1)
#
#   nsim = 1000
#   pvals <- numeric(nsim)
#   score_stats <- numeric(nsim)
#   for(iter in 1:nsim){
#     print(iter)
#     Y <- matrix(NA,ncol = J, nrow = n)
#
#     # set.seed(iter)
#     for(i in 1:n){
#       for(j in 1:J){
#         temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#         # Y[i,j] <- rpois(1, lambda = temp_mean)*rbinom(1,1,0.4)
#         Y[i,j] <- rnbinom(1, mu = temp_mean,size = 1.5)*rbinom(1,1,0.4)
#       }
#     }
#
#     B <- suppressMessages(
#       emuFit_micro_penalized(X,Y)
#      )
#
#     Y <- B$Y_augmented
#     B <- B$B
#
#     # set.seed(iter)
#     score_result <- #suppressMessages(
#       try(micro_score_test(Y= Y,
#                            X = X,
#                            B= B,
#                            null_k = 2,
#                            null_j = J/2,
#                            tolerance = 0.1,
#                            rho_scaling = 10,
#                            rho_init = 1,
#                            gap_tolerance = 1e-5,
#                            constr_init_tol = 1e6,
#                            constr_tol_scaling = 0.8,
#                            maxit =5000,
#                            constraint_fn = function(x){
#                              pseudohuber_center(x,0.1)},
#                            constraint_grad_fn = function(x){
#                              dpseudohuber_center_dx(x,0.1)
#                            },
#                            constraint_hess_fn = function(x,ind_1,ind_2){
#                              hess_pseudohuber_center(x,0.1,ind_1,ind_2)
#                            }))
#     # )
#
#     pvals[iter] <- score_result$pval
#     # print(pvals[iter])
#
#     score_stats[iter] <- score_result$score_stat
#     if(iter %%10 ==0){
#       # hist(score_stats[1:iter],xlim = c(-10,10),freq= FALSE,breaks= 50)
#       # x = seq(-10,10,by = .01)
#       # lines(x,dt(x,n - 1),lty= 2)
#       hist(log(score_stats[1:iter]),freq= FALSE,breaks= 50,
#            xlim = c(-20,5))
#       abline(v = log(qchisq(0.95,1)),lty = 2, col = "red")
#       x <- exp(seq(-20,5,0.1))
#       lines(log(x),x*dchisq(x,1),lty = 3)
#       print(signif(mean(pvals[1:iter]<=0.05)))
#     }
#
#   }
#
#
#    })
# # #
# # test_that("Test score test with Poisson data under null, J > n", {
# #   set.seed(4323)
# #   n <- 8
# #   J <- 20
# #   X <- cbind(1,rep(c(0,1),each = n/2))
# #   z <- rnorm(n) + 8
# #   b0 <- rnorm(J)
# #   b1 <- seq(-1.5,1.5,length.out = J)
# #   b1[J/2 + 0:1] <- 0
# #   b <- rbind(b0,b1)
# #
# #   nsim = 1000
# #   pvals <- numeric(nsim)
# #   score_stats <- numeric(nsim)
# #   for(iter in 1:nsim){
# #     print(iter)
# #     Y <- matrix(NA,ncol = J, nrow = n)
# #
# #     for(i in 1:n){
# #       for(j in 1:J){
# #         temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
# #         # Y[i,j] <- rpois(1, lambda = temp_mean)
# #         Y[i,j] <- rnbinom(1, mu = temp_mean,size = 1.5)
# #       }
# #     }
# #
# #     B <- suppressMessages(emuFit_micro(X,Y))
# #     # set.seed(iter)
# #     score_result <- #suppressMessages(
# #       try(micro_score_test(Y= Y,
# #                                      X = X,
# #                                      B= B,
# #                                      null_k = 2,
# #                                      null_j = J/2,
# #                                      tolerance = 0.1,
# #                                      rho_scaling = 10,
# #                                      rho_init = 10,
# #                                      gap_tolerance = 1e-4,
# #                            constr_init_tol = 1e5,
# #                            maxit =5000,
# #                                      constraint_fn = function(x){
# #                                        pseudohuber_center(x,0.1)},
# #                                      constraint_grad_fn = function(x){
# #                                        dpseudohuber_center_dx(x,0.1)
# #                                      },
# #                                      constraint_hess_fn = function(x,ind_1,ind_2){
# #                                        hess_pseudohuber_center(x,0.1,ind_1,ind_2)
# #                                      }))
# #
# #     # diagnostic <- lagrangian_fit(X,
# #     #                            Y,
# #     #                            B = B,   maxit = 1000,
# #     #                constraint_fn = function(x){
# #     #                  pseudohuber_center(x,0.1)},
# #     #                constraint_grad_fn = function(x){
# #     #                  dpseudohuber_center_dx(x,0.1)
# #     #                },
# #     #                constraint_hess_fn = function(x,ind_1,ind_2){
# #     #                  hess_pseudohuber_center(x,0.1,ind_1,ind_2)
# #     #                },
# #     # j_constr = 10,
# #     # k_constr = 2,
# #     #                            max_step = 0.5,
# #     #                            rho_init = .1,
# #     #                            rho_scaling = 10,
# #     #                            tolerance = 1e-1,
# #     #                            gap_tolerance = 1e-4,
# #     #                            step_ratio = 0.5,
# #     #                            c1 = 1e-4,
# #     #                            verbose = TRUE,
# #     #                            constr_init_tol = 1e5,
# #     #                            constr_tol_scaling = 0.8,
# #     #                            n_to_update = 5)
# #       # )
# # # diagnostic %>%
# # #   group_by(k,j) %>%
# # #   mutate(value = value - value[iter == max(iter)]) %>%
# # #   filter(iter >90) %>%
# # #   ggplot() +
# # #   geom_line(aes(x = iter,y = value, group = interaction(k,j))) +
# # #   geom_point(aes(x = iter, y = value, color = stepsize))
# # #
# # # diagnostic %>%
# # #   ggplot() +
# # #   geom_line(aes(x = iter, y = stepsize, group = interaction(k,j))) +
# # #   scale_y_log10()
# #
# #     pvals[iter] <- score_result$pval
# #     # print(pvals[iter])
# #
# #     score_stats[iter] <- score_result$score_stat
# #     if(iter %%10 ==0){
# #     # hist(score_stats[1:iter],xlim = c(-10,10),freq= FALSE,breaks= 50)
# #     # x = seq(-10,10,by = .01)
# #     # lines(x,dt(x,n - 1),lty= 2)
# #       hist(score_stats[1:iter],xlim = c(0,10),freq= FALSE,breaks= 50)
# #     x <- seq(0,10,by = 0.01)
# #     lines(x,dchisq(x,1))
# #     print(signif(mean(pvals[1:iter]<=0.05)))
# # }
# #
# #   }
# #
# #   nb_8_20 <- list("score_stats" = score_stats,
# #                   "pvals" = pvals)
# #
# #
# #   set.seed(4323)
# #   n <- 8
# #   J <- 20
# #   X <- cbind(1,rep(c(0,1),each = n/2))
# #   z <- rnorm(n) + 8
# #   b0 <- rnorm(J)
# #   b1 <- seq(-1.5,1.5,length.out = J)
# #   b1[J/2 + 0:1] <- 0
# #   b <- rbind(b0,b1)
# #
# #   nsim = 1000
# #   pvals <- numeric(nsim)
# #   score_stats <- numeric(nsim)
# #   for(iter in 1:nsim){
# #     print(iter)
# #     Y <- matrix(NA,ncol = J, nrow = n)
# #
# #     for(i in 1:n){
# #       for(j in 1:J){
# #         temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
# #         Y[i,j] <- rpois(1, lambda = temp_mean)
# #         # Y[i,j] <- rnbinom(1, mu = temp_mean,size = 1.5)
# #       }
# #     }
# #
# #     B <- suppressMessages(emuFit_micro(X,Y))
# #     # set.seed(iter)
# #     score_result <- #suppressMessages(
# #       try(micro_score_test(Y= Y,
# #                            X = X,
# #                            B= B,
# #                            null_k = 2,
# #                            null_j = J/2,
# #                            tolerance = 1e-3,
# #                            rho_scaling = 10,
# #                            rho_init = 10,
# #                            gap_tolerance = 1e-4,
# #                            constr_init_tol = 1e5,
# #                            maxit =1000,
# #                            constraint_fn = function(x){
# #                              pseudohuber_center(x,0.1)},
# #                            constraint_grad_fn = function(x){
# #                              dpseudohuber_center_dx(x,0.1)
# #                            },
# #                            constraint_hess_fn = function(x,ind_1,ind_2){
# #                              hess_pseudohuber_center(x,0.1,ind_1,ind_2)
# #                            }))
# #
# #     # diagnostic <- lagrangian_fit(X,
# #     #                            Y,
# #     #                            B = B,   maxit = 1000,
# #     #                constraint_fn = function(x){
# #     #                  pseudohuber_center(x,0.1)},
# #     #                constraint_grad_fn = function(x){
# #     #                  dpseudohuber_center_dx(x,0.1)
# #     #                },
# #     #                constraint_hess_fn = function(x,ind_1,ind_2){
# #     #                  hess_pseudohuber_center(x,0.1,ind_1,ind_2)
# #     #                },
# #     # j_constr = 10,
# #     # k_constr = 2,
# #     #                            max_step = 0.5,
# #     #                            rho_init = .1,
# #     #                            rho_scaling = 10,
# #     #                            tolerance = 1e-1,
# #     #                            gap_tolerance = 1e-4,
# #     #                            step_ratio = 0.5,
# #     #                            c1 = 1e-4,
# #     #                            verbose = TRUE,
# #     #                            constr_init_tol = 1e5,
# #     #                            constr_tol_scaling = 0.8,
# #     #                            n_to_update = 5)
# #     # )
# #     # diagnostic %>%
# #     #   group_by(k,j) %>%
# #     #   mutate(value = value - value[iter == max(iter)]) %>%
# #     #   filter(iter >90) %>%
# #     #   ggplot() +
# #     #   geom_line(aes(x = iter,y = value, group = interaction(k,j))) +
# #     #   geom_point(aes(x = iter, y = value, color = stepsize))
# #     #
# #     # diagnostic %>%
# #     #   ggplot() +
# #     #   geom_line(aes(x = iter, y = stepsize, group = interaction(k,j))) +
# #     #   scale_y_log10()
# #
# #     pvals[iter] <- score_result$pval
# #     # print(pvals[iter])
# #
# #     score_stats[iter] <- score_result$score_stat
# #     if(iter %%10 ==0){
# #       # hist(score_stats[1:iter],xlim = c(-10,10),freq= FALSE,breaks= 50)
# #       # x = seq(-10,10,by = .01)
# #       # lines(x,dt(x,n - 1),lty= 2)
# #       hist(score_stats[1:iter],xlim = c(0,10),freq= FALSE,breaks= 50)
# #       x <- seq(0,10,by = 0.01)
# #       lines(x,dchisq(x,1))
# #       print(signif(mean(pvals[1:iter]<=0.05)))
# #     }
# #
# #   }
# #  pois_8_20 <- list("score_stats" = score_stats,
# #                   "pvals" = pvals)
# #
# #
# #
# #  set.seed(4323)
# #  n <- 1000
# #  J <- 20
# #  X <- cbind(1,rep(c(0,1),each = n/2))
# #  z <- rnorm(n) + 8
# #  b0 <- rnorm(J)
# #  b1 <- seq(-1.5,1.5,length.out = J)
# #  b1[J/2 + 0:1] <- 0
# #  b <- rbind(b0,b1)
# #
# #  nsim = 1000
# #  pvals <- numeric(nsim)
# #  score_stats <- numeric(nsim)
# #  for(iter in 1:nsim){
# #    print(iter)
# #    Y <- matrix(NA,ncol = J, nrow = n)
# #
# #    for(i in 1:n){
# #      for(j in 1:J){
# #        temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
# #        Y[i,j] <- rpois(1, lambda = temp_mean)
# #        # Y[i,j] <- rnbinom(1, mu = temp_mean,size = 1.5)
# #      }
# #    }
# #
# #    B <- suppressMessages(emuFit_micro(X,Y))
# #    # set.seed(iter)
# #    score_result <- #suppressMessages(
# #      try(micro_score_test(Y= Y,
# #                           X = X,
# #                           B= B,
# #                           null_k = 2,
# #                           null_j = J/2,
# #                           tolerance = 1,
# #                           rho_scaling = 10,
# #                           rho_init = 10,
# #                           gap_tolerance = 1e-4,
# #                           constr_init_tol = 1e5,
# #                           maxit =1000,
# #                           constraint_fn = function(x){
# #                             pseudohuber_center(x,0.1)},
# #                           constraint_grad_fn = function(x){
# #                             dpseudohuber_center_dx(x,0.1)
# #                           },
# #                           constraint_hess_fn = function(x,ind_1,ind_2){
# #                             hess_pseudohuber_center(x,0.1,ind_1,ind_2)
# #                           }))
# #
# #    # diagnostic <- lagrangian_fit(X,
# #    #                            Y,
# #    #                            B = B,   maxit = 1000,
# #    #                constraint_fn = function(x){
# #    #                  pseudohuber_center(x,0.1)},
# #    #                constraint_grad_fn = function(x){
# #    #                  dpseudohuber_center_dx(x,0.1)
# #    #                },
# #    #                constraint_hess_fn = function(x,ind_1,ind_2){
# #    #                  hess_pseudohuber_center(x,0.1,ind_1,ind_2)
# #    #                },
# #    # j_constr = 10,
# #    # k_constr = 2,
# #    #                            max_step = 0.5,
# #    #                            rho_init = .1,
# #    #                            rho_scaling = 10,
# #    #                            tolerance = 1e-1,
# #    #                            gap_tolerance = 1e-4,
# #    #                            step_ratio = 0.5,
# #    #                            c1 = 1e-4,
# #    #                            verbose = TRUE,
# #    #                            constr_init_tol = 1e5,
# #    #                            constr_tol_scaling = 0.8,
# #    #                            n_to_update = 5)
# #    # )
# #    # diagnostic %>%
# #    #   group_by(k,j) %>%
# #    #   mutate(value = value - value[iter == max(iter)]) %>%
# #    #   filter(iter >90) %>%
# #    #   ggplot() +
# #    #   geom_line(aes(x = iter,y = value, group = interaction(k,j))) +
# #    #   geom_point(aes(x = iter, y = value, color = stepsize))
# #    #
# #    # diagnostic %>%
# #    #   ggplot() +
# #    #   geom_line(aes(x = iter, y = stepsize, group = interaction(k,j))) +
# #    #   scale_y_log10()
# #
# #    pvals[iter] <- score_result$pval
# #    # print(pvals[iter])
# #
# #    score_stats[iter] <- score_result$score_stat
# #    if(iter %%10 ==0){
# #      # hist(score_stats[1:iter],xlim = c(-10,10),freq= FALSE,breaks= 50)
# #      # x = seq(-10,10,by = .01)
# #      # lines(x,dt(x,n - 1),lty= 2)
# #      hist(score_stats[1:iter],xlim = c(0,10),freq= FALSE,breaks= 50)
# #      x <- seq(0,10,by = 0.01)
# #      lines(x,dchisq(x,1))
# #      print(signif(mean(pvals[1:iter]<=0.05)))
# #    }
# #
# #  }
# #  pois_500_20 <- list("score_stats" = score_stats,
# #                    "pvals" = pvals)
# # #
# # hist(pvals,freq = FALSE,breaks = 100,xlim = c(0,1))
# #  mean(pvals <=0.05)
# #
# # hist(score_stats,breaks =200,freq=FALSE)
# # x = seq(-10,10,by = .01)
# # lines(x,dt(x,n - 1),lty= 2,col="red")
# #
# # qqplot(score_stats,rt(1e5,n-1))
# # abline(a= 0, b=1, col  ="red")
# # sd(score_stats)
# #
# #
# # })
#
# # test_that("Test score test with Poisson data under null, J > n, huber constraint", {
# #   set.seed(4323)
# #   n <- 10
# #   J <- 20
# #   X <- cbind(1,rep(c(-.5,.5),each = n/2))
# #   z <- rnorm(n) +3
# #   b0 <- rnorm(J)
# #   b1 <- seq(-1.5,1.5,length.out = J)
# #   b1[J/2 + 0:1] <- 0
# #   b <- rbind(b0,b1)
# #
# #   nsim = 1e4
# #   pvals <- numeric(nsim)
# #   score_stats <- numeric(nsim)
# #   for(iter in 1:nsim){
# #     print(iter)
# #     Y <- matrix(NA,ncol = J, nrow = n)
# #
# #     for(i in 1:n){
# #       for(j in 1:J){
# #         temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
# #         # Y[i,j] <- rpois(1, lambda = temp_mean)
# #         Y[i,j] <- rnbinom(1, mu = temp_mean,size = 1.5)
# #       }
# #     }
# #
# #     score_result <- try(micro_score_test(Y= Y,
# #                                          X = X,
# #                                          B= NULL,
# #                                          constraint_fn = huber_center,
# #                                          null_k = 2,
# #                                          null_j = J/2,
# #                                          tolerance = 1e-4,
# #                                          constraint_type = "huber",
# #                                          huber_param = 1))
# #
# #
# #     pvals[iter] <- score_result$pval
# #     # print(pvals[iter])
# #
# #     score_stats[iter] <- score_result$score_stat
# #     if(iter %%10 ==0){
# #       # hist(score_stats[1:iter],xlim = c(-10,10),freq= FALSE,breaks= 50)
# #       # x = seq(-10,10,by = .01)
# #       # lines(x,dt(x,n - 1),lty= 2)
# #       hist(score_stats[1:iter],xlim = c(0,10),freq= FALSE,breaks= 50)
# #       x <- seq(0,10,by = 0.01)
# #       lines(x,dchisq(x,1))
# #       print(signif(mean(pvals[1:iter]<=0.05)))
# #     }
# #
# #   }
# #
# #   hist(pvals,freq = FALSE,breaks = 100)
# #   mean(pvals <=0.05)
# #
# #   hist(score_stats,breaks =200,freq=FALSE)
# #   x = seq(-10,10,by = .01)
# #   lines(x,dt(x,n - 1),lty= 2,col="red")
# #
# #   qqplot(score_stats,rt(1e5,n-1))
# #   abline(a= 0, b=1, col  ="red")
# #   sd(score_stats)
# #
# #
# # })
# #
# # test_that("Test score test with Poisson data under null, larger n", {
# #   set.seed(4323)
# #   n <- 20
# #   J <- 80
# #   X <- cbind(1,rep(c(-.5,.5),each = n/2))
# #   z <- rnorm(n) +3
# #   b0 <- rnorm(J)
# #   b1 <- seq(-1.5,1.5,length.out = J)
# #   b1[J/2 + 0:1] <- 0
# #   b <- rbind(b0,b1)
# #
# #   nsim = 1e4
# #   pvals <- numeric(nsim)
# #   score_stats <- numeric(nsim)
# #   for(iter in 1:nsim){
# #     print(iter)
# #     Y <- matrix(NA,ncol = J, nrow = n)
# #
# #     for(i in 1:n){
# #       for(j in 1:J){
# #         temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
# #         Y[i,j] <- rpois(1, lambda = temp_mean)
# #         # Y[i,j] <- rnbinom(1, mu = temp_mean,size = 1.5)
# #       }
# #     }
# #
# #     score_result <- try(micro_score_test(Y= Y,
# #                                          X = X,
# #                                          B= NULL,
# #                                          constraint_fn = mean,
# #                                          null_k = 2,
# #                                          null_j = J/2,
# #                                          tolerance = 1e-5,
# #                                          type = "correct_score",
# #                                          cov_est_method = "plugin",
# #                                          t_se_method= "from_cov",
# #                                          loo = FALSE))
# #
# #
# #     pvals[iter] <- score_result$pval
# #     # print(pvals[iter])
# #
# #     score_stats[iter] <- score_result$score_stat
# #     if(iter %%10 ==0){
# #       # hist(score_stats[1:iter],xlim = c(-10,10),freq= FALSE,breaks= 50)
# #       # x = seq(-10,10,by = .01)
# #       # lines(x,dt(x,n - 1),lty= 2)
# #       hist(score_stats[1:iter],xlim = c(0,10),freq= FALSE,breaks= 50)
# #       x <- seq(0,10,by = 0.01)
# #       lines(x,dchisq(x,1))
# #       print(signif(mean(pvals[1:iter]<=0.05)))
# #     }
# #
# #   }
# #
# # hist(pvals,freq = FALSE,breaks = 100)
# #   mean(pvals <=0.05)
# #
# #   hist(score_stats,breaks =200,freq=FALSE)
# #   x = seq(-10,10,by = .01)
# #   lines(x,dt(x,n - 1),lty= 2,col="red")
# #
# #   qqplot(score_stats,rt(1e5,n-1))
# #   abline(a= 0, b=1, col  ="red")
# #   sd(score_stats)
# #
# #
# # })
# #
# # test_that("Test score test with Poisson data under null, J > n", {
# #   set.seed(4323)
# #   n <- 6
# #   J <- 20
# #   X <- cbind(1,rep(c(-.5,.5),each = n/2))
# #   z <- rnorm(n) +3
# #   b0 <- rnorm(J)
# #   b1 <- seq(-1.5,1.5,length.out = J)
# #   b1[J/2 + 0:1] <- 0
# #   b <- rbind(b0,b1)
# #
# #   results <- vector(4,mode = "list")
# #
# #   nsim = 1e4
# #   pvals <- numeric(nsim)
# #   score_stats <- numeric(nsim)
# #   for(iter in 1:nsim){
# #     print(iter)
# #     Y <- matrix(NA,ncol = J, nrow = n)
# #
# #     for(i in 1:n){
# #       for(j in 1:J){
# #         temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
# #         # Y[i,j] <- rpois(1, lambda = temp_mean)
# #         Y[i,j] <- rnbinom(1, mu = temp_mean,size = 2)
# #       }
# #     }
# #
# #     score_result <- try(micro_score_test(Y= Y,
# #                                          X = X,
# #                                          B= NULL,
# #                                          constraint_fn = mean,
# #                                          null_k = 2,
# #                                          null_j = J/2,
# #                                          tolerance = 1e-5,
# #                                          maxit = 50,
# #                                          type = "t",
# #                                          cov_est_method = "crossproduct",
# #                                          t_se_method= "from_cov",
# #                                          loo = FALSE))
# #
# #
# #     pvals[iter] <- score_result$pval
# #     # print(pvals[iter])
# #
# #     score_stats[iter] <- score_result$score_stat
# #     if(iter %%10 ==0){
# #       # hist(score_stats[1:iter],xlim = c(-10,10),freq= FALSE,breaks= 50)
# #       # x = seq(-10,10,by = .01)
# #       # lines(x,dt(x,n - 1),lty= 2)
# #       hist(score_stats[1:iter],xlim = c(0,10),freq= FALSE,breaks= 50)
# #       x <- seq(0,10,by = 0.01)
# #       lines(x,dchisq(x,1),col = "red")
# #       print(signif(mean(pvals[1:iter]<=0.05)))
# #     }
# #
# #   }
# #
# #   results[[1]] <- list("score_stats" = score_stats,
# #                        "pvals" = pvals)
# #
# #   pvals <- numeric(nsim)
# #   score_stats <- numeric(nsim)
# #   for(iter in 1:nsim){
# #     print(iter)
# #     Y <- matrix(NA,ncol = J, nrow = n)
# #
# #     for(i in 1:n){
# #       for(j in 1:J){
# #         temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
# #         # Y[i,j] <- rpois(1, lambda = temp_mean)
# #         Y[i,j] <- rnbinom(1, mu = temp_mean,size = 1.5)
# #       }
# #     }
# #
# #     score_result <- try(micro_score_test(Y= Y,
# #                                          X = X,
# #                                          B= NULL,
# #                                          constraint_fn = mean,
# #                                          null_k = 2,
# #                                          null_j = J/2,
# #                                          tolerance = 1e-5,
# #                                          type = "t",
# #                                          cov_est_method = "crossproduct",
# #                                          t_se_method= "from_cov",
# #                                          loo = FALSE))
# #
# #
# #     pvals[iter] <- score_result$pval
# #     # print(pvals[iter])
# #
# #     score_stats[iter] <- score_result$score_stat
# #     if(iter %%10 ==0){
# #       # hist(score_stats[1:iter],xlim = c(-10,10),freq= FALSE,breaks= 50)
# #       # x = seq(-10,10,by = .01)
# #       # lines(x,dt(x,n - 1),lty= 2)
# #       hist(score_stats[1:iter],xlim = c(0,10),freq= FALSE,breaks= 50)
# #       x <- seq(0,10,by = 0.01)
# #       lines(x,dchisq(x,1))
# #       print(signif(mean(pvals[1:iter]<=0.05)))
# #     }
# #
# #   }
# #
# #   results[[2]] <- list("score_stats" = score_stats,
# #                        "pvals" = pvals)
# #   set.seed(4323)
# #   n <- 10
# #   J <- 80
# #   X <- cbind(1,rep(c(-.5,.5),each = n/2))
# #   z <- rnorm(n) +3
# #   b0 <- rnorm(J)
# #   b1 <- seq(-1.5,1.5,length.out = J)
# #   b1[J/2 + 0:1] <- 0
# #   b <- rbind(b0,b1)
# #   pvals <- numeric(nsim)
# #   score_stats <- numeric(nsim)
# #   for(iter in 1:nsim){
# #     print(iter)
# #     Y <- matrix(NA,ncol = J, nrow = n)
# #
# #     for(i in 1:n){
# #       for(j in 1:J){
# #         temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
# #         Y[i,j] <- rpois(1, lambda = temp_mean)
# #         # Y[i,j] <- rnbinom(1, mu = temp_mean,size = 1.5)
# #       }
# #     }
# #
# #     score_result <- try(micro_score_test(Y= Y,
# #                                          X = X,
# #                                          B= NULL,
# #                                          constraint_fn = mean,
# #                                          null_k = 2,
# #                                          null_j = J/2,
# #                                          tolerance = 1e-5,
# #                                          type = "t",
# #                                          cov_est_method = "crossproduct",
# #                                          t_se_method= "from_cov",
# #                                          loo = FALSE))
# #
# #
# #     pvals[iter] <- score_result$pval
# #     # print(pvals[iter])
# #
# #     score_stats[iter] <- score_result$score_stat
# #     if(iter %%10 ==0){
# #       # hist(score_stats[1:iter],xlim = c(-10,10),freq= FALSE,breaks= 50)
# #       # x = seq(-10,10,by = .01)
# #       # lines(x,dt(x,n - 1),lty= 2)
# #       hist(score_stats[1:iter],xlim = c(0,10),freq= FALSE,breaks= 50)
# #       x <- seq(0,10,by = 0.01)
# #       lines(x,dchisq(x,1))
# #       print(signif(mean(pvals[1:iter]<=0.05)))
# #     }
# #
# #   }
# #
# #   results[[3]] <- list("score_stats" = score_stats,"pvals" = pvals)
# #
# #   pvals <- numeric(nsim)
# #   score_stats <- numeric(nsim)
# #   for(iter in 1:nsim){
# #     print(iter)
# #     Y <- matrix(NA,ncol = J, nrow = n)
# #
# #     for(i in 1:n){
# #       for(j in 1:J){
# #         temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
# #         # Y[i,j] <- rpois(1, lambda = temp_mean)
# #         Y[i,j] <- rnbinom(1, mu = temp_mean,size = 1.5)
# #       }
# #     }
# #
# #     score_result <- try(micro_score_test(Y= Y,
# #                                          X = X,
# #                                          B= NULL,
# #                                          constraint_fn = mean,
# #                                          null_k = 2,
# #                                          null_j = J/2,
# #                                          tolerance = 1e-5,
# #                                          type = "t",
# #                                          cov_est_method = "crossproduct",
# #                                          t_se_method= "from_cov",
# #                                          loo = FALSE))
# #
# #
# #     pvals[iter] <- score_result$pval
# #     # print(pvals[iter])
# #
# #     score_stats[iter] <- score_result$score_stat
# #     if(iter %%10 ==0){
# #       # hist(score_stats[1:iter],xlim = c(-10,10),freq= FALSE,breaks= 50)
# #       # x = seq(-10,10,by = .01)
# #       # lines(x,dt(x,n - 1),lty= 2)
# #       hist(score_stats[1:iter],xlim = c(0,10),freq= FALSE,breaks= 50)
# #       x <- seq(0,10,by = 0.01)
# #       lines(x,dchisq(x,1))
# #       print(signif(mean(pvals[1:iter]<=0.05)))
# #     }
# #
# #   }
# #
# #   results[[4]] <- list("score_stats" = score_stats,"pvals" = pvals)
# #
# # })
# #
# # data.frame("J" = rep(c(20,80),each= 2),
# #            "n" = 10,
# #            "Distribution" = rep(c("Poisson","NB with size = 1.5")),
# #            "empirical_alpha"= sapply(1:4,function(i) mean(results[[i]]$pval<=0.05))) %>%
# #   ggplot() +
# #   geom_point(aes(x = J,y = empirical_alpha,color = Distribution)) +
# #   ylim(c(0,.1)) +
# #   geom_abline(aes(intercept = 0.05,slope= 0),linetype = 2)+
# #   xlab("Number of Taxa") +
# #   ylab("Empirical Rejection Rate at 0.05 Level") +
# #   theme_bw()
# #
# # rbind(
# #   data.frame("Distribution" = "Poisson",
# #              "pvals" = results[[1]]$pvals,
# #              "J" = 20),
# #   data.frame("Distribution" = "NB with size = 1.5",
# #              "pvals" = results[[2]]$pvals,
# #              "J" = 20),
# #   data.frame("Distribution" = "Poisson",
# #              "pvals" = results[[3]]$pvals,
# #              "J" = 80),
# #   data.frame("Distribution" = "NB with size = 1.5",
# #              "pvals" = results[[4]]$pvals,
# #              "J" = 80)
# # ) %>%
# #   ggplot(aes(sample = pvals)) +
# #   stat_qq(distribution = qunif,size = 0.5,alpha = 0.01) +
# #   geom_abline(aes(intercept = 0, slope = 1),linetype = 2, color= "red") +
# #   facet_grid(J~Distribution) +
# #   theme_bw()
# #
# #
# # test_that("Test score test with Poisson data under null", {
# #   set.seed(4323)
# #   n <- 20
# #   J <- 10
# #   X <- cbind(1,rep(c(-.5,.5),each = n/2))
# #   z <- rnorm(n) +1
# #   b0 <- rnorm(J)
# #   b1 <- seq(-1.5,1.5,length.out = J)
# #   b1[J/2 + 0:1] <- 0
# #   b <- rbind(b0,b1)
# #
# #   nsim = 1e4
# #   pvals <- numeric(nsim)
# #   score_stats <- numeric(nsim)
# #   for(iter in 1:nsim){
# #     print(iter)
# #     Y <- matrix(NA,ncol = J, nrow = n)
# #
# #     for(i in 1:n){
# #       for(j in 1:J){
# #         temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
# #         Y[i,j] <- rpois(1, lambda = temp_mean)
# #         # Y[i,j] <- rnbinom(1, mu = temp_mean,size = 1.5)
# #       }
# #     }
# #
# #     pf <- suppressMessages(emuFit_micro_penalized(X[rowSums(Y)>0,],
# #                   Y[rowSums(Y)>0,],
# #                   B = NULL,
# #                   constraint_fn = function(x) x[J],
# #                   maxit = 500,
# #                   tolerance = 1e-3,
# #                   collect_iterations = FALSE))
# #
# #     score_result <- try(micro_score_test(Y= pf$Y_augmented,
# #                                          X = X[rowSums(Y)>0,],
# #                                          B= NULL,
# #                                          constraint_fn = mean,
# #                                          null_k = 2,
# #                                          null_j = J/2,
# #                                          tolerance = 1e-5,
# #                                          type = "correct_score",
# #                                          cov_est_method = "plugin",
# #                                          t_se_method= "from_cov",
# #                                          loo = FALSE))
# #
# #
# #     pvals[iter] <- score_result$pval
# #     # print(pvals[iter])
# #
# #     score_stats[iter] <- score_result$score_stat
# #     if(iter %%10 ==0){
# #       # hist(score_stats[1:iter],xlim = c(-10,10),freq= FALSE,breaks= 50)
# #       # x = seq(-10,10,by = .01)
# #       # lines(x,dt(x,n - 1),lty= 2)
# #       hist(score_stats[1:iter],xlim = c(0,10),freq= FALSE,breaks= 50)
# #       x <- seq(0,10,by = 0.01)
# #       lines(x,dchisq(x,1))
# #       print(signif(mean(pvals[1:iter]<=0.05)))
# #     }
# #
# #   }
# #
# #     hist(pvals,freq = FALSE,breaks = 100)
# #   mean(pvals <=0.05)
# #
# #   hist(score_stats,breaks =200,freq=FALSE)
# #   x = seq(-10,10,by = .01)
# #   lines(x,dt(x,n - 1),lty= 2,col="red")
# #
# #   qqplot(score_stats,rt(1e5,n-1))
# #   abline(a= 0, b=1, col  ="red")
# #   sd(score_stats)
# #
# #
# # })
# #
# #
