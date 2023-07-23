#
#
# emuFit_micro_constrained <-
#   function(X,
#            Y,
#            B = NULL,
#            constraint_fn = NULL,
#            constraint_grad_fn = NULL,
#            j_constr,
#            k_constr,
#            maxit = 250,
#            tolerance = 1e-1,
#            info_reg = 0.5,
#            max_step = 0.5,
#            c1 = 1e-4,
#            verbose = TRUE
#   ){
#     n <- nrow(Y)
#     J <- ncol(Y)
#     p <- ncol(X)
#
#     if(is.null(B)){
#       B <- matrix(0,
#                   nrow = p,
#                   ncol = J)
#     }
#
#     if(is.null(constraint_fn)){
#       constraint_fn <- pseudohuber_center
#       constraint_grad_fn <- dpseudohuber_center_dx
#     }
#
#     #directly enforce constraint at starting value of B
#     B[k_constr,j_constr] <- constraint_fn(B[k_constr,-j_constr])
#
#     for(k in 1:p){
#       B[k,] <- B[k,] - constraint_fn(B[k,])
#     }
#
#     z <- update_z_no_wts(Y = Y,
#                          X = X,
#                          B = B)
#
#     converged <- FALSE
#     iter <- 1
#     lls <- numeric(0)
#
#     # B_list <- vector(0,mode = "list")
#     B_cup <- B_cup_from_B(B)
#
#     B_cup_df <- data.frame(matrix(ncol= nrow(B_cup) , nrow = maxit +1))
#     B_cup_df$iter <- 0:maxit
#     B_cup_df[1,1:(p*J)] <- as.numeric(as.matrix(B_cup))
#
#
#     while(!converged){
#
#       for(j in 1:J){
#         update <- micro_fisher_constr(X = X,
#                                       Y = Y,
#                                       B = B,
#                                       z = z,
#                                       j = j,
#                                       j_constr = j_constr,
#                                       k_constr = k_constr,
#                                       constraint_fn = constraint_fn,
#                                       constraint_grad_fn = constraint_grad_fn,
#                                       stepsize = max_step,
#                                       c1 = c1)
#
#
#         B[,j] <- B[,j] + update
#
#         if(j != j_constr){
#           B[k_constr,j_constr] <- constraint_fn(B[k_constr,-j_constr])
#         }
#
#         for(k in 1:p){
#           B[k,] <- B[k,] - constraint_fn(B[k,])
#         }
#
#
#         z <- update_z_no_wts(Y = Y,X = X,B = B)
#         logmu <- X%*%B + matrix(z, ncol = 1)%*%matrix(1,ncol = J, nrow = 1)
#         ll <- round(sum(Y*logmu - exp(logmu)),2)
#         lls <- c(lls,ll)
#
#       }
#
#       # for(k in 1:p){
#       #   B[k,] <- B[k,] - constraint_fn(B[k,])
#       # }
#       # #
#       # #
#       # z <- update_z_no_wts(Y = Y,X = X,B = B)
#
#       B_cup_df[iter + 1,1:(p*J)] <- B_cup_from_B(B)
#
#       # B_list <- append(B_list,list(B))
#
#       lagrangian_gr <-  constr_lagrangian(Y = Y,
#                                           X = X,
#                                           B = B,
#                                           constraint_fn = constraint_fn,
#                                           constraint_grad_fn = constraint_grad_fn,
#                                           k_constr = k_constr,
#                                           j_constr = j_constr)
#
#       deriv_norm <- sum(sqrt(lagrangian_gr^2))
#       nll <- length(lls)
#       n_to_plot <- min(nll,100)
#       if(lls[nll] - lls[nll - 1]<0){
#         warning("ll decreased!")
#       }
#
#       plot(1:n_to_plot + nll - n_to_plot,lls[(nll - n_to_plot + 1):nll],type = "l")
#
#
#       if(verbose){
#         message("Norm of lagrangian derivative ",signif(deriv_norm,4) )
#       }
#
#       if(deriv_norm < tolerance){
#         converged <- TRUE
#       }
#
#       iter <- iter + 1
#
#       if(iter > maxit){
#         converged <- TRUE
#         if(verbose){
#           message("Exiting optimization; maximum number of iterations reached.")
#         }
#       }
#
#     }
#
#
#
#     return(B)
#
#   }
#
# #
# #
# # emuFit_micro_constrained <-
# #   function(X,
# #            Y,
# #            B = NULL,
# #            constraint_fn = NULL,
# #            constraint_grad_fn = NULL,
# #            j_constr,
# #            k_constr,
# #            maxit = 250,
# #            tolerance = 1e-1,
# #            info_reg = 0.5,
# #            max_step = 0.5,
# #            # step_decay = TRUE,
# #            abs_max_step = 0.5,
# #            c1 = 1e-4,
# #            step_ratio = 0.5,
# #            verbose = TRUE
# #            ){
# #     n <- nrow(Y)
# #     J <- ncol(Y)
# #     p <- ncol(X)
# #     orig_step <- max_step
# #
# #     if(is.null(B)){
# #       B <- matrix(0,
# #                   nrow = p,
# #                   ncol = J)
# #     }
# #
# #     if(is.null(constraint_fn)){
# #       constraint_fn <- pseudohuber_center
# #       constraint_grad_fn <- dpseudohuber_center_dx
# #     }
# #
# #     #directly enforce constraint at starting value of B
# #
# #     B[k_constr,] <- constraint_fn(B[k_constr,-j_constr])
# #     B[k_constr,j_constr] <- constraint_fn(B[k_constr,-j_constr])
# #
# #     for(k in 1:p){
# #       B[k,] <- B[k,] - constraint_fn(B[k,])
# #     }
# #
# #     z <- update_z_no_wts(Y = Y,
# #                          X = X,
# #                          B = B)
# #
# #     converged <- FALSE
# #     iter <- 1
# #     lls <- lgs <-  numeric(0)
# #
# #     # B_list <- vector(0,mode = "list")
# #     B_cup <- B_cup_from_B(B)
# #
# #     B_cup_df <- data.frame(matrix(ncol= nrow(B_cup) , nrow = maxit +1))
# #     B_cup_df$iter <- 0:maxit
# #     B_cup_df[1,1:(p*J)] <- as.numeric(as.matrix(B_cup))
# #
# #
# #     while(!converged){
# #
# #       for(j in 1:J){
# #         update <- micro_fisher_constr(X = X,
# #                                       Y = Y,
# #                                       B = B,
# #                                       z = z,
# #                                       j = j,
# #                                       j_constr = j_constr,
# #                                       k_constr = k_constr,
# #                                       constraint_fn = constraint_fn,
# #                                       constraint_grad_fn = constraint_grad_fn,
# #                                       stepsize = max_step,
# #                                       step_ratio = step_ratio,
# #                                       c1 = c1)
# #
# #       B[,j] <- B[,j] + update
# #
# #
# #       if(j != j_constr){
# #         B[k_constr,j_constr] <- constraint_fn(B[k_constr,-j_constr])
# #       }
# #
# #         for(k in 1:p){
# #           B[k,] <- B[k,] - constraint_fn(B[k,])
# #         }
# #
# #
# #         z <- update_z_no_wts(Y = Y,X = X,B = B)
# #         logmu <- X%*%B + matrix(z, ncol = 1)%*%matrix(1,ncol = J, nrow = 1)
# #         ll <- round(sum(Y*logmu - exp(logmu)),2)
# #         lls <- c(lls,ll)
# #
# #         lg <-   constr_lagrangian(Y = Y,
# #                                                        X = X,
# #                                                        B = B,
# #                                                        constraint_fn = constraint_fn,
# #                                                        constraint_grad_fn = constraint_grad_fn,
# #                                                        k_constr = k_constr,
# #                                                        j_constr = j_constr)
# #
# #         lgs <- c(lgs,sqrt(sum(lg^2)))
# #
# #       }
# #
# #
# #       # if(step_decay){
# #       #   max_step <- orig_step/(1 + 0.01*iter)
# #       # }
# #       # for(k in 1:p){
# #       #   B[k,] <- B[k,] - constraint_fn(B[k,])
# #       # }
# #       # #
# #       # #
# #       # z <- update_z_no_wts(Y = Y,X = X,B = B)
# #
# #       B_cup_df[iter + 1,1:(p*J)] <- B_cup_from_B(B)
# #
# #       # B_list <- append(B_list,list(B))
# #
# #      lagrangian_gr <-  constr_lagrangian(Y = Y,
# #                                     X = X,
# #                                     B = B,
# #                                     constraint_fn = constraint_fn,
# #                                     constraint_grad_fn = constraint_grad_fn,
# #                                     k_constr = k_constr,
# #                                     j_constr = j_constr)
# #
# #      deriv_norm <- sqrt(sum(lagrangian_gr^2))
# #      nlg <- length(lgs)
# #      n_to_plot <- min(nlg,100)
# #      if(lgs[nlg] - lgs[nlg - 1]>0){
# #        warning("norm of gradient increased!")
# #      }
# #
# #      plot(1:n_to_plot + nlg - n_to_plot,lgs[(nlg - n_to_plot + 1):nlg],type = "l")
# #
# #
# #        if(verbose){
# #          message("Norm of lagrangian derivative ",signif(deriv_norm,4) )
# #        }
# #
# #      if(deriv_norm < tolerance){
# #        converged <- TRUE
# #      }
# #
# #      iter <- iter + 1
# #
# #      if(iter > maxit){
# #        converged <- TRUE
# #        if(verbose){
# #        message("Exiting optimization; maximum number of iterations reached.")
# #          }
# #      }
# #
# #     }
# #
# #
# #
# #     return(B)
# #
# #   }
# #
# # this_iter <- iter
# # # #
# # B_names <- matrix("a",nrow = nrow(B),ncol = ncol(B))
# # for(j in 1:J){
# #   for(k in 1:p){
# #     B_names[k,j] <- paste("k",k,"j",j,collapse = "_",sep = "_")
# #   }
# # }
# #
# # J <- ncol(B)
# # p <- nrow(B)
# # B_name_cup <- matrix("",nrow = J*p,ncol = 1)
# # for(j in 1:J){
# #   B_name_cup[(j - 1)*p + 1:p] <- B_names[,j]
# # }
# # colnames(B_cup_df)[1:(p*J)] <- as.character(B_name_cup)
# # B_cup_df %>%
# #   filter(iter<5) %>%
# #   filter(iter<this_iter) %>%
# #   # filter(iter>100) %>%
# #   pivot_longer(-iter) %>%
# #   group_by(name)%>%
# #   # mutate(value = value - value[iter == (this_iter - 1)],) %>%
# #   ggplot() +
# #   geom_line(aes(x = iter,y = value, group = name,color = name)) +
# #   theme_bw()
# # #
# # #
# # #
# # # #
# # #
# # #
# # # B_cup_df %>%
# # #   filter(iter<314) %>%
# # #   filter(iter>250) %>%
# # #   pivot_longer(-iter) %>%
# # #   group_by(name)%>%
# # #   summarize(sd_trail = sd(value)) %>%
# # #   with(cbind(name,sd_trail)[order(sd_trail,decreasing = TRUE),])
# #
# #
# #
# #
# # # n <- nrow(Y)
# # # J <- ncol(Y)
# # # p <- ncol(X)
# # #
# # #
# # # j_constr <- constrained_kj[2]
# # # k_constr <- constrained_kj[1]
# # #
# # # if(is.null(B)){
# # #   B <- matrix(0,
# # #               nrow = p,
# # #               ncol = J)
# # # }
# # #
# # # z <- update_z_no_wts(Y,X,B)
# # #
# # # log_mean <- X%*%B +
# # #   matrix(z,ncol = 1)%*%matrix(1,ncol = J, nrow = 1)
# # #
# # # lls <- sum(Y*log_mean - exp(log_mean))
# # # granular_lls <- lls
# # # steps <- "start"
# # #
# # # converged <- FALSE
# # #
# # # iter <- 1
# # # plot(as.numeric(B),ylim = c(-10,10))
# # # while(!converged){
# # #   points(as.numeric(B))
# # #   B_old <- B
# # #   # print(iter)
# # #   if(iter %% 50 ==0){
# # #     # print(B)
# # #   }
# # #   for(j in 1:J){
# # #
# # #     dg_dBj <- get_dg_dBj(B = B,
# # #                          j = j,
# # #                          constraint_type = constraint_type,
# # #                          k_star =  k_constr,
# # #                          huber_param = huber_param)
# # # #
# # #     update = micro_fisher_constr(X = X,
# # #                                  Y = Y,
# # #                                  B = B,
# # #                                  z = z,
# # #                                  j = j,
# # #                                  jstar = j_constr,
# # #                                  kstar = k_constr,
# # #                                  dg_dBj = dg_dBj)
# # #
# # #     B[,j] <- B[,j] + update
# # #     if(j != j_constr){
# # #       if(constraint_type != "huber"){
# # #         B[k_constr,j_constr] <- constraint_fn(B[k_constr,-j_constr])
# # #       } else{
# # #         B[k_constr,j_constr] <- huber_center(B[k_constr,-j_constr],
# # #                                              d = huber_param)
# # #       }
# # #     }
# # #
# # #     z <- update_z_no_wts(Y,X,B)
# # #
# # #   }
# # #   for(k in 1:p){
# # #     if(constraint_type != "huber"){
# # #       B[k,] <- B[k,] - constraint_fn(B[k,]) #constraint_fn(B[k_constr,-j_constr])
# # #     } else{
# # #       B[k,] <- B[k,] - huber_center(B[k,],d = huber_param)
# # #     }
# # #   }
# # #   print(max(abs(B_old - B)))
# # #   if(max(abs(B_old - B)) < tolerance){
# # #     converged <- TRUE
# # #   }
# # #   iter <- iter + 1
# # #   print(iter)
# # #   if(iter> maxit){
# # #     converge <- TRUE
# # #   }
# # # }
# # # return(B)
# # #
# # # return(list(B = B,
# # #             ll = lls[length(lls)],
# # #             lls = lls))
# #
#
# #
# # emuFit_micro_constrained <-
# #   function(X,
# #            Y,
# #            B = NULL,
# #            constraint_fn = NULL,
# #            constraint_grad_fn = NULL,
# #            j_constr,
# #            k_constr,
# #            maxit = 250,
# #            tolerance = 1e-1,
# #            info_reg = 0.5,
# #            max_step = 0.5,
# #            c1 = 1e-4,
# #            verbose = TRUE
# #   ){
# #     n <- nrow(Y)
# #     J <- ncol(Y)
# #     p <- ncol(X)
# #
# #     if(is.null(B)){
# #       B <- matrix(0,
# #                   nrow = p,
# #                   ncol = J)
# #     }
# #
# #     if(is.null(constraint_fn)){
# #       constraint_fn <- pseudohuber_center
# #       constraint_grad_fn <- dpseudohuber_center_dx
# #     }
# #
# #     #directly enforce constraint at starting value of B
# #     B[k_constr,j_constr] <- constraint_fn(B[k_constr,-j_constr])
# #
# #     for(k in 1:p){
# #       B[k,] <- B[k,] - constraint_fn(B[k,])
# #     }
# #
# #     z <- update_z_no_wts(Y = Y,
# #                          X = X,
# #                          B = B)
# #
# #     converged <- FALSE
# #     iter <- 1
# #     lls <- numeric(0)
# #
# #     # B_list <- vector(0,mode = "list")
# #     B_cup <- B_cup_from_B(B)
# #
# #     B_cup_df <- data.frame(matrix(ncol= nrow(B_cup) , nrow = maxit +1))
# #     B_cup_df$iter <- 0:maxit
# #     B_cup_df[1,1:(p*J)] <- as.numeric(as.matrix(B_cup))
# #
# #     lagrangian_gr <-  constr_lagrangian(Y = Y,
# #                                         X = X,
# #                                         B = B,
# #                                         constraint_fn = constraint_fn,
# #                                         constraint_grad_fn = constraint_grad_fn,
# #                                         k_constr = k_constr,
# #                                         j_constr = j_constr)
# #
# #     deriv_norm <- sqrt(sum(lagrangian_gr^2))
# #     while(!converged){
# #
# #       for(j in 1:J){
# #
# #         local_grad <- Inf
# #         mini_iter <- 1
# #       while(sqrt(sum(local_grad^2))>1 & mini_iter<25){
# #
# #         update <- micro_fisher_constr(X = X,
# #                                       Y = Y,
# #                                       B = B,
# #                                       z = z,
# #                                       j = j,
# #                                       j_constr = j_constr,
# #                                       k_constr = k_constr,
# #                                       constraint_fn = constraint_fn,
# #                                       constraint_grad_fn = constraint_grad_fn,
# #                                       stepsize = max_step,
# #                                       c1 = c1)
# #
# #         if(mini_iter==1){
# #           old_update <- update
# #         } else{
# #           arch_update <- update
# #           update <- update + old_update
# #           old_update <- arch_update
# #         }
# #
# #
# #
# #         B[,j] <- B[,j] + update
# #
# #
# #
# #
# #         # print(update)
# #         #
# #         mini_iter <- mini_iter + 1
# #         # print(mini_iter)
# #
# #
# #
# #         if(j != j_constr){
# #           B[k_constr,j_constr] <- constraint_fn(B[k_constr,-j_constr])
# #         }
# #
# #         for(k in 1:p){
# #           B[k,] <- B[k,] - constraint_fn(B[k,])
# #         }
# #
# #
# #         z <- update_z_no_wts(Y = Y,X = X,B = B)
# #
# #
# #         if(j != j_constr){
# #           local_grad <- dljaug_dBj(Y,
# #                                    X,
# #                                    B,
# #                                    z,
# #                                    constraint_fn,
# #                                    constraint_grad_fn,
# #                                    j = j,
# #                                    j_constr,
# #                                    k_constr
# #           )
# #
# #           # ll_fn <- function(Bj){
# #           #   temp_B <- B
# #           #   temp_B[,j] <- Bj
# #           #   temp_B[,j_constr] <- constraint_fn(temp_B[k_constr,-j_constr])
# #           #   temp_means <- exp(X%*%temp_B + matrix(z,ncol = 1)%*%matrix(1,nrow = 1, ncol = J))
# #           #   return(sum(Y*log(temp_means) - temp_means))
# #           # }
# #           #
# #           # numDeriv::grad(ll_fn,B[,j])
# #
# #
# #         } else{
# #           means <- exp(X%*%B[,j,drop = FALSE] + z)
# #           local_grad <- t(X[,-k_constr,drop = FALSE])%*%(Y[,j, drop = FALSE] - matrix(means,ncol = 1))
# #         }
# #         # print(sqrt(sum(local_grad^2)))
# #       }
# #
# #
# #
# #
# #         logmu <- X%*%B + matrix(z, ncol = 1)%*%matrix(1,ncol = J, nrow = 1)
# #         ll <- round(sum(Y*logmu - exp(logmu)),2)
# #         lls <- c(lls,ll)
# #
# #       }
# #
# #       # for(k in 1:p){
# #       #   B[k,] <- B[k,] - constraint_fn(B[k,])
# #       # }
# #       # #
# #       # #
# #       # z <- update_z_no_wts(Y = Y,X = X,B = B)
# #
# #       B_cup_df[iter + 1,1:(p*J)] <- B_cup_from_B(B)
# #
# #       # B_list <- append(B_list,list(B))
# #
# #       lagrangian_gr <-  constr_lagrangian(Y = Y,
# #                                           X = X,
# #                                           B = B,
# #                                           constraint_fn = constraint_fn,
# #                                           constraint_grad_fn = constraint_grad_fn,
# #                                           k_constr = k_constr,
# #                                           j_constr = j_constr)
# #
# #       deriv_norm <- sqrt(sum(lagrangian_gr^2))
# #       # nll <- length(lls)
# #       # n_to_plot <- min(nll,100)
# #       # if(lls[nll] - lls[nll - 1]<0){
# #       #   warning("ll decreased!")
# #       # }
# #
# #       # plot(1:n_to_plot + nll - n_to_plot,lls[(nll - n_to_plot + 1):nll],type = "l")
# #
# #
# #       if(verbose){
# #         message("Norm of lagrangian derivative ",signif(deriv_norm,4) )
# #       }
# #
# #       if(deriv_norm < tolerance){
# #         converged <- TRUE
# #       }
# #
# #       iter <- iter + 1
# #
# #       if(iter > maxit){
# #         converged <- TRUE
# #         if(verbose){
# #           message("Exiting optimization; maximum number of iterations reached.")
# #         }
# #       }
# #
# #     }
# #
# #
# #
# #     return(B)
# #
# #   }
#
# #
# #
# # emuFit_micro_constrained <-
# #   function(X,
# #            Y,
# #            B = NULL,
# #            constraint_fn = NULL,
# #            constraint_grad_fn = NULL,
# #            j_constr,
# #            k_constr,
# #            maxit = 250,
# #            tolerance = 1e-1,
# #            info_reg = 0.5,
# #            max_step = 0.5,
# #            # step_decay = TRUE,
# #            abs_max_step = 0.5,
# #            c1 = 1e-4,
# #            step_ratio = 0.5,
# #            verbose = TRUE
# #            ){
# #     n <- nrow(Y)
# #     J <- ncol(Y)
# #     p <- ncol(X)
# #     orig_step <- max_step
# #
# #     if(is.null(B)){
# #       B <- matrix(0,
# #                   nrow = p,
# #                   ncol = J)
# #     }
# #
# #     if(is.null(constraint_fn)){
# #       constraint_fn <- pseudohuber_center
# #       constraint_grad_fn <- dpseudohuber_center_dx
# #     }
# #
# #     #directly enforce constraint at starting value of B
# #
# #     B[k_constr,] <- constraint_fn(B[k_constr,-j_constr])
# #     B[k_constr,j_constr] <- constraint_fn(B[k_constr,-j_constr])
# #
# #     for(k in 1:p){
# #       B[k,] <- B[k,] - constraint_fn(B[k,])
# #     }
# #
# #     z <- update_z_no_wts(Y = Y,
# #                          X = X,
# #                          B = B)
# #
# #     converged <- FALSE
# #     iter <- 1
# #     lls <- lgs <-  numeric(0)
# #
# #     # B_list <- vector(0,mode = "list")
# #     B_cup <- B_cup_from_B(B)
# #
# #     B_cup_df <- data.frame(matrix(ncol= nrow(B_cup) , nrow = maxit +1))
# #     B_cup_df$iter <- 0:maxit
# #     B_cup_df[1,1:(p*J)] <- as.numeric(as.matrix(B_cup))
# #
# #
# #     while(!converged){
# #
# #       for(j in 1:J){
# #         update <- micro_fisher_constr(X = X,
# #                                       Y = Y,
# #                                       B = B,
# #                                       z = z,
# #                                       j = j,
# #                                       j_constr = j_constr,
# #                                       k_constr = k_constr,
# #                                       constraint_fn = constraint_fn,
# #                                       constraint_grad_fn = constraint_grad_fn,
# #                                       stepsize = max_step,
# #                                       step_ratio = step_ratio,
# #                                       c1 = c1)
# #
# #       B[,j] <- B[,j] + update
# #
# #
# #       if(j != j_constr){
# #         B[k_constr,j_constr] <- constraint_fn(B[k_constr,-j_constr])
# #       }
# #
# #         for(k in 1:p){
# #           B[k,] <- B[k,] - constraint_fn(B[k,])
# #         }
# #
# #
# #         z <- update_z_no_wts(Y = Y,X = X,B = B)
# #         logmu <- X%*%B + matrix(z, ncol = 1)%*%matrix(1,ncol = J, nrow = 1)
# #         ll <- round(sum(Y*logmu - exp(logmu)),2)
# #         lls <- c(lls,ll)
# #
# #         lg <-   constr_lagrangian(Y = Y,
# #                                                        X = X,
# #                                                        B = B,
# #                                                        constraint_fn = constraint_fn,
# #                                                        constraint_grad_fn = constraint_grad_fn,
# #                                                        k_constr = k_constr,
# #                                                        j_constr = j_constr)
# #
# #         lgs <- c(lgs,sqrt(sum(lg^2)))
# #
# #       }
# #
# #
# #       # if(step_decay){
# #       #   max_step <- orig_step/(1 + 0.01*iter)
# #       # }
# #       # for(k in 1:p){
# #       #   B[k,] <- B[k,] - constraint_fn(B[k,])
# #       # }
# #       # #
# #       # #
# #       # z <- update_z_no_wts(Y = Y,X = X,B = B)
# #
# #       B_cup_df[iter + 1,1:(p*J)] <- B_cup_from_B(B)
# #
# #       # B_list <- append(B_list,list(B))
# #
# #      lagrangian_gr <-  constr_lagrangian(Y = Y,
# #                                     X = X,
# #                                     B = B,
# #                                     constraint_fn = constraint_fn,
# #                                     constraint_grad_fn = constraint_grad_fn,
# #                                     k_constr = k_constr,
# #                                     j_constr = j_constr)
# #
# #      deriv_norm <- sqrt(sum(lagrangian_gr^2))
# #      nlg <- length(lgs)
# #      n_to_plot <- min(nlg,100)
# #      if(lgs[nlg] - lgs[nlg - 1]>0){
# #        warning("norm of gradient increased!")
# #      }
# #
# #      plot(1:n_to_plot + nlg - n_to_plot,lgs[(nlg - n_to_plot + 1):nlg],type = "l")
# #
# #
# #        if(verbose){
# #          message("Norm of lagrangian derivative ",signif(deriv_norm,4) )
# #        }
# #
# #      if(deriv_norm < tolerance){
# #        converged <- TRUE
# #      }
# #
# #      iter <- iter + 1
# #
# #      if(iter > maxit){
# #        converged <- TRUE
# #        if(verbose){
# #        message("Exiting optimization; maximum number of iterations reached.")
# #          }
# #      }
# #
# #     }
# #
# #
# #
# #     return(B)
# #
# #   }
# #
# this_iter <- iter
# # #
# B_names <- matrix("a",nrow = nrow(B),ncol = ncol(B))
# for(j in 1:J){
#   for(k in 1:p){
#     B_names[k,j] <- paste("k",k,"j",j,collapse = "_",sep = "_")
#   }
# }
#
# J <- ncol(B)
# p <- nrow(B)
# B_name_cup <- matrix("",nrow = J*p,ncol = 1)
# for(j in 1:J){
#   B_name_cup[(j - 1)*p + 1:p] <- B_names[,j]
# }
# colnames(B_cup_df)[1:(p*J)] <- as.character(B_name_cup)
# B_cup_df %>%
#   # filter(iter<5) %>%
#   filter(iter<this_iter) %>%
#   # filter(iter>100) %>%
#   pivot_longer(-iter) %>%
#   group_by(name)%>%
#   # mutate(value = value - value[iter == (this_iter - 1)],) %>%
#   ggplot() +
#   geom_line(aes(x = iter,y = value, group = name,color = name)) +
#   theme_bw()
# # #
# # #
# # #
# # # #
# # #
# # #
# # # B_cup_df %>%
# # #   filter(iter<314) %>%
# # #   filter(iter>250) %>%
# # #   pivot_longer(-iter) %>%
# # #   group_by(name)%>%
# # #   summarize(sd_trail = sd(value)) %>%
# # #   with(cbind(name,sd_trail)[order(sd_trail,decreasing = TRUE),])
# #
# #
# #
# #
# # # n <- nrow(Y)
# # # J <- ncol(Y)
# # # p <- ncol(X)
# # #
# # #
# # # j_constr <- constrained_kj[2]
# # # k_constr <- constrained_kj[1]
# # #
# # # if(is.null(B)){
# # #   B <- matrix(0,
# # #               nrow = p,
# # #               ncol = J)
# # # }
# # #
# # # z <- update_z_no_wts(Y,X,B)
# # #
# # # log_mean <- X%*%B +
# # #   matrix(z,ncol = 1)%*%matrix(1,ncol = J, nrow = 1)
# # #
# # # lls <- sum(Y*log_mean - exp(log_mean))
# # # granular_lls <- lls
# # # steps <- "start"
# # #
# # # converged <- FALSE
# # #
# # # iter <- 1
# # # plot(as.numeric(B),ylim = c(-10,10))
# # # while(!converged){
# # #   points(as.numeric(B))
# # #   B_old <- B
# # #   # print(iter)
# # #   if(iter %% 50 ==0){
# # #     # print(B)
# # #   }
# # #   for(j in 1:J){
# # #
# # #     dg_dBj <- get_dg_dBj(B = B,
# # #                          j = j,
# # #                          constraint_type = constraint_type,
# # #                          k_star =  k_constr,
# # #                          huber_param = huber_param)
# # # #
# # #     update = micro_fisher_constr(X = X,
# # #                                  Y = Y,
# # #                                  B = B,
# # #                                  z = z,
# # #                                  j = j,
# # #                                  jstar = j_constr,
# # #                                  kstar = k_constr,
# # #                                  dg_dBj = dg_dBj)
# # #
# # #     B[,j] <- B[,j] + update
# # #     if(j != j_constr){
# # #       if(constraint_type != "huber"){
# # #         B[k_constr,j_constr] <- constraint_fn(B[k_constr,-j_constr])
# # #       } else{
# # #         B[k_constr,j_constr] <- huber_center(B[k_constr,-j_constr],
# # #                                              d = huber_param)
# # #       }
# # #     }
# # #
# # #     z <- update_z_no_wts(Y,X,B)
# # #
# # #   }
# # #   for(k in 1:p){
# # #     if(constraint_type != "huber"){
# # #       B[k,] <- B[k,] - constraint_fn(B[k,]) #constraint_fn(B[k_constr,-j_constr])
# # #     } else{
# # #       B[k,] <- B[k,] - huber_center(B[k,],d = huber_param)
# # #     }
# # #   }
# # #   print(max(abs(B_old - B)))
# # #   if(max(abs(B_old - B)) < tolerance){
# # #     converged <- TRUE
# # #   }
# # #   iter <- iter + 1
# # #   print(iter)
# # #   if(iter> maxit){
# # #     converge <- TRUE
# # #   }
# # # }
# # # return(B)
# # #
# # # return(list(B = B,
# # #             ll = lls[length(lls)],
# # #             lls = lls))
# #
