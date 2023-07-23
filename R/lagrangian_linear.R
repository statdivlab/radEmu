# lagrangian_fit <- function(X,
#                            Y,
#                            B = NULL,
#                            constraint_fn = NULL,
#                            constraint_grad_fn = NULL,
#                            j_constr,
#                            k_constr,
#                            maxit = 250,
#                            rho_init = 10,
#                            rho_scaling = 1.5,
#                            tolerance = 1e-1,
#                            gap_tolerance = 1e-6,
#                            info_reg = 0.5,
#                            max_step = 0.5,
#                            c1 = 1e-4,
#                            verbose = TRUE){
# J <- ncol(Y)
# p <- ncol(X)
# n <- nrow(Y)
# z <- update_z_no_wts(Y,X,B)
#
#
# gap <- B[k_constr,j_constr] - constraint_fn(B[k_constr,])
#
# rho <- rho_init
#
# u <- rho*gap
#
# # ek_constr <- matrix(0,nrow = p, ncol = 1)
# # ek_constr[k_constr,] <- 1
# # ek_constr_out <- ek_constr%*%t(ek_constr)
#
#
#
# B0 <- B
#
# # message("David! lin_const() should be its own function with tests and stuff!")
# # lin_const <- function(x,con_gr,x0){
# #   return(
# #     as.numeric(
# #       constraint_fn(B0[k_constr,]) +
# #    crossprod(con_gr,matrix(x - x0,ncol = 1))
# #     )
# #   )
# # }
#
# feas_gap <- B[k_constr,j_constr] - constraint_fn(B[k_constr,])
#
# B_cups <- data.frame(value = B_cup[,1],
#                      j = rep(1:J,each = 2),
#                      k = rep(1:2, J),
#                      rho = NA,
#                      u = NA,
#                      iter = 0)
# iter <- 1
# while(abs(feas_gap)>1e-8){
#
# constraint_grad_at_B0 <- constraint_grad_fn(B[k_constr,])
# constraint_fn_at_B0 <- constraint_fn(B[k_constr,])
# lgns <- numeric(0)
# lag_gr_norm <- Inf
# B0 <- B
# while(lag_gr_norm>0.1){
# for(j in 1:J){
# #   means <- exp(X%*%B[,j,drop = FALSE] + z)
# #   message("David! lag_fn() should *definitely* be its own function with
# # tests and it should *definitely* not rely on variables not passed to it!")
# #   lag_fn <- function(Bj){
# #     temp_B <- B
# #     temp_B[,j] <- Bj
# #     temp_gap <- (temp_B[k_constr,j_constr] -
# #                    lin_const(temp_B[k_constr,], con_gr, B0[k_constr,]))
# #     log_means <- X%*%temp_B + matrix(z, ncol = 1)%*%matrix(1, nrow = 1, ncol = J)
# #     return(-sum(Y*log_means - exp(log_means)) +
# #       u*temp_gap +
# #       rho*temp_gap^2)
# #   }
# # message("Same thing for lag_gr() -- you know this.")
# #   lag_gr <- -t(X)%*%(Y[,j] - means) +
# #     u*(ifelse(j == j_constr,1,0) - con_gr[j])*ek_constr +
# #     rho*(B[k_constr,j_constr] -
# #            lin_const(B[k_constr,],
# #                      con_gr,
# #                      old_B[k_constr,]))*(ifelse(j == j_constr,1,0) - con_gr[j])*ek_constr
# # message("And lag_h() -- own function + tests!")
# #   lag_h <- t(X)%*%diag(as.numeric(means))%*%X +
# #     rho*ek_constr_out*(ifelse(j == j_constr,1,0) - con_gr[j])^2
#
#   curr_vals <-
#     linearized_aug_lag_z(X = X,
#                        Y = Y,
#                        j = j,
#                        j_constr = j_constr,
#                        k_constr = k_constr,
#                        B = B,
#                        z = z,
#                        B0 = B0,
#                        u = u,
#                        rho = rho,
#                        constraint_fn_at_B0 = constraint_fn_at_B0,
#                        constraint_grad_at_B0 = constraint_grad_at_B0,
#                        constraint_fn = constraint_fn,
#                        constraint_grad_fn = constraint_grad_fn,
#                        compute_gradient = TRUE,
#                        compute_hessian = TRUE)
#
#   # update <- try(qr.solve(curr_vals$hess,curr_vals$gr),
#   #               silent = TRUE)
#   # if(inherits(update,"try-error")){
#     update <- try(qr.solve(curr_vals$hess+ diag(rep(sqrt(mean(curr_vals$gr^2)),p)),curr_vals$gr))
#   # }
#
#   update_norm <- max(abs(update))
#   # if(update_norm>0.1){
#   #   update <- 0.1*update/update_norm
#   # }
#   stepsize <- 0.5
#   accept <- FALSE
#
#
#   armijo <- c1*sum(lag_gr*update)
#   while(!accept){
#     prop_B <- B
#     prop_B[,j] <- B[,j] - stepsize*update
#     prop_val <-    linearized_aug_lag_z(X = X,
#                                         Y = Y,
#                                         j = j,
#                                         j_constr = j_constr,
#                                         k_constr = k_constr,
#                                         B = prop_B,
#                                         z = z,
#                                         B0 = B0,
#                                         u = u,
#                                         rho = rho,
#                                         constraint_fn_at_B0 = constraint_fn_at_B0,
#                                         constraint_grad_at_B0 = constraint_grad_at_B0,
#                                         constraint_fn = constraint_fn,
#                                         constraint_grad_fn = constraint_grad_fn)
#
#     accept <- prop_val$value <= curr_vals$value+ armijo*stepsize
#     stepsize <- stepsize/2
#   }
#   B[,j] <-   prop_B[,j]
#   for(k in 1:p){
#     B[k,] <- B[k,] - constraint_fn(B[k,])
#   }
#   z <- update_z_no_wts(Y,X,B)
#
# }
#
#   B_cups <- rbind(B_cups,
#                   data.frame(value = B_cup_from_B(B)[,1],
#                        j = rep(1:J,each = 2),
#                        k = rep(1:2, J),
#                        rho = rho,
#                        u = u,
#                        iter = iter))
#
#   iter <- iter + 1
#
#
#
# # message("see how terrible this is?")
#   lag_gr <- do.call(c,
#                     lapply(1:J,
#                            function(j)
#                              linearized_aug_lag_z(X = X,
#                                                   Y = Y,
#                                                   j = j,
#                                                   j_constr = j_constr,
#                                                   k_constr = k_constr,
#                                                   B = prop_B,
#                                                   z = z,
#                                                   B0 = B0,
#                                                   u = u,
#                                                   rho = rho,
#                                                   constraint_fn_at_B0 = constraint_fn_at_B0,
#                                                   constraint_grad_at_B0 = constraint_grad_at_B0,
#                                                   constraint_fn = constraint_fn,
#                                                   constraint_grad_fn = constraint_grad_fn,
#                                                   compute_gradient = TRUE)$gr
#
#                            )
#   )
#
#
#     # do.call(c,
#     #
#     #         # lapply(1:J,
#     #         #   function(j){
#     #         #     means <- exp(X%*%B[,j,drop = FALSE] + z);
#     #         #   return(-t(X)%*%(Y[,j] - means) +
#     #         #   u*(ifelse(j == j_constr,1,0) - con_gr[j])*ek_constr +
#     #         #   rho*(B[k_constr,j_constr] -
#     #         #   lin_const(B[k_constr,],
#     #         #   con_gr,
#     #         #   B0[k_constr,]))*
#     #         #     (ifelse(j == j_constr,1,0) - con_gr[j])*ek_constr)
#     #         #     }
#     #         #   )
#     # )
#
#   lag_gr_norm <- sqrt(sum(lag_gr^2))
#   lgns <- c(lgns,lag_gr_norm)
#
#   # print(lag_gr_norm)
#
# }
#
# # step <- B - B0
# # step_norm <- max(abs(step))
# #
# # if(step_norm>0.1){
# #   step <- 0.1*step/step_norm
# # }
# # B <- B0 + step
#
# feas_gap <- B[k_constr,j_constr] - constraint_fn(B[k_constr,])
# print(feas_gap)
# u <- u + rho*feas_gap
# us <- c(us,u)
# rho <- rho*2
#
#
#
# }
#
# }
#
# }
# B_cups %>%
#   ggplot() +
#   # geom_point(aes(x = iter, y= asinh(value),color = u))
#   geom_line(aes(x = iter,y = asinh(value), group = interaction(k,j)))
