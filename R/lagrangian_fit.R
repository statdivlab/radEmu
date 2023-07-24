
lagrangian_fit <- function(X,
                           Y,
                           B = NULL,
                           constraint_fn = NULL,
                           constraint_grad_fn = NULL,
                           j_constr,
                           k_constr,
                           maxit = 250,
                           rho_init = 10,
                           rho_scaling = 2,
                           tolerance = 1e-1,
                           gap_tolerance = 1e-6,
                           c1 = 1e-4,
                           verbose = TRUE){
  J <- ncol(Y)
  p <- ncol(X)
  n <- nrow(Y)
  z <- update_z_no_wts(Y,X,B)


  gap <- B[k_constr,j_constr] - constraint_fn(B[k_constr,])

  rho <- rho_init

  u <- rho*gap
  # us <- u

  # ek_constr <- matrix(0,nrow = p, ncol = 1)
  # ek_constr[k_constr,] <- 1
  # ek_constr_out <- ek_constr%*%t(ek_constr)



  B0 <- B

  # message("David! lin_const() should be its own function with tests and stuff!")
  # lin_const <- function(x,con_gr,x0){
  #   return(
  #     as.numeric(
  #       constraint_fn(B0[k_constr,]) +
  #    crossprod(con_gr,matrix(x - x0,ncol = 1))
  #     )
  #   )
  # }

  feas_gap <- B[k_constr,j_constr] - constraint_fn(B[k_constr,])
# 
  # B_cups <- data.frame(value = B_cup_from_B(B)[,1],
  #                      j = rep(1:J,each = 2),
  #                      k = rep(1:2, J),
  #                      rho = NA,
  #                      u = NA,
  #                      iter = 0)
  
  
  lag_gr <- do.call(c,
                    lapply(1:J,
                           function(j)
                             linearized_aug_lag_z(X = X,
                                                  Y = Y,
                                                  j = j,
                                                  j_constr = j_constr,
                                                  k_constr = k_constr,
                                                  B = B,
                                                  z = z,
                                                  B0 = B0,
                                                  u = u,
                                                  rho = rho,
                                                  constraint_fn = constraint_fn,
                                                  constraint_grad_fn = constraint_grad_fn,
                                                  compute_gradient = TRUE)$gr
                           
                    )
  )
  iter <- 1
  while(abs(feas_gap)>gap_tolerance){

    constraint_grad_at_B0 <- constraint_grad_fn(B[k_constr,])
    constraint_fn_at_B0 <- constraint_fn(B[k_constr,])
    lgns <- numeric(0)
    lag_gr_norm <- Inf
    B0 <- B
    while(lag_gr_norm>tolerance){
      for(j in 1:J){
        #   means <- exp(X%*%B[,j,drop = FALSE] + z)
        #   message("David! lag_fn() should *definitely* be its own function with
        # tests and it should *definitely* not rely on variables not passed to it!")
        #   lag_fn <- function(Bj){
        #     temp_B <- B
        #     temp_B[,j] <- Bj
        #     temp_gap <- (temp_B[k_constr,j_constr] -
        #                    lin_const(temp_B[k_constr,], con_gr, B0[k_constr,]))
        #     log_means <- X%*%temp_B + matrix(z, ncol = 1)%*%matrix(1, nrow = 1, ncol = J)
        #     return(-sum(Y*log_means - exp(log_means)) +
        #       u*temp_gap +
        #       rho*temp_gap^2)
        #   }
        # message("Same thing for lag_gr() -- you know this.")
        #   lag_gr <- -t(X)%*%(Y[,j] - means) +
        #     u*(ifelse(j == j_constr,1,0) - con_gr[j])*ek_constr +
        #     rho*(B[k_constr,j_constr] -
        #            lin_const(B[k_constr,],
        #                      con_gr,
        #                      old_B[k_constr,]))*(ifelse(j == j_constr,1,0) - con_gr[j])*ek_constr
        # message("And lag_h() -- own function + tests!")
        #   lag_h <- t(X)%*%diag(as.numeric(means))%*%X +
        #     rho*ek_constr_out*(ifelse(j == j_constr,1,0) - con_gr[j])^2

        curr_vals <-
          linearized_aug_lag_z(X = X,
                               Y = Y,
                               j = j,
                               j_constr = j_constr,
                               k_constr = k_constr,
                               B = B,
                               z = z,
                               B0 = B0,
                               u = u,
                               rho = rho,
                               constraint_fn = constraint_fn,
                               constraint_grad_fn = constraint_grad_fn,
                               compute_gradient = TRUE,
                               compute_hessian = TRUE)
        
        # message("test derivatives being returned by linearized_aug_lag_z")
# 
#         al_func <- function(Bj){
#           temp_B <- B
#           temp_B[,j] <- Bj
#           linearized_aug_lag_z(X = X,
#                                      Y = Y,
#                                      j = j,
#                                      j_constr = j_constr,
#                                      k_constr = k_constr,
#                                      B = temp_B,
#                                      z = z,
#                                      B0 = B0,
#                                      u = u,
#                                      rho = rho,
#                                      constraint_fn = constraint_fn,
#                                      constraint_grad_fn = constraint_grad_fn,
#                                      compute_gradient = FALSE,
#                                      compute_hessian = FALSE)[[1]]
#         }
#         #
#         al_func(B[,j])
#         #
#         nd <- numDeriv::grad(al_func,B[,j])
#         nh <- numDeriv::hessian(al_func,B[,j])
#         ad <- curr_vals$gr
#         ah <- curr_vals$hess
# 
#         ad/nd
        
     
        
        
        
        # log_means <- X%*%B + matrix(z,ncol = 1)%*%matrix(1,nrow = 1, ncol = J)
        # 
        # ll <- sum(Y*log_means - exp(log_means))
        # 
        # gap <- B[null_k,null_j] - constraint_fn(B[null_k,])
        # 
        # -ll + u*gap + 0.5*rho*gap^2

        update <- try(qr.solve(curr_vals$hess,curr_vals$gr),
                      silent = TRUE)
        reg <- 0.1
        while(inherits(update,"try-error")){
        update <- try(qr.solve(curr_vals$hess+ diag(reg*rep(sqrt(mean(curr_vals$gr^2)),p)),curr_vals$gr))
        reg <- 10*reg
        }

        update_norm <- max(abs(update))
        # if(update_norm>0.1){
        #   update <- 0.1*update/update_norm
        # }
        stepsize <- 0.5
        accept <- FALSE

        # message("finding stepsize for category ", j)
        armijo <- c1*sum(curr_vals$gr*update)
        while(!accept){
          prop_B <- B
          prop_B[,j] <- B[,j] - stepsize*update
          prop_val <-    linearized_aug_lag_z(X = X,
                                              Y = Y,
                                              j = j,
                                              j_constr = j_constr,
                                              k_constr = k_constr,
                                              B = prop_B,
                                              z = z,
                                              B0 = B0,
                                              u = u,
                                              rho = rho,
                                              constraint_fn = constraint_fn,
                                              constraint_grad_fn = constraint_grad_fn)

          accept <- prop_val$value <= curr_vals$value+ armijo*stepsize
          stepsize <- stepsize*0.9
        }
        B[,j] <-   prop_B[,j]
        for(k in 1:p){
          B[k,] <- B[k,] - constraint_fn(B[k,])
        }
        z <- update_z_no_wts(Y,X,B)

      }
      # message("step accepted")

      # B_cups <- rbind(B_cups,
      #                 data.frame(value = B_cup_from_B(B)[,1],
      #                            j = rep(1:J,each = 2),
      #                            k = rep(1:2, J),
      #                            rho = rho,
      #                            u = u,
      #                            iter = iter))

      iter <- iter + 1



      # message("see how terrible this is?")
      lag_gr <- do.call(c,
                        lapply(1:J,
                               function(j)
                                 linearized_aug_lag_z(X = X,
                                                      Y = Y,
                                                      j = j,
                                                      j_constr = j_constr,
                                                      k_constr = k_constr,
                                                      B = B,
                                                      z = z,
                                                      B0 = B0,
                                                      u = u,
                                                      rho = rho,
                                                      constraint_fn = constraint_fn,
                                                      constraint_grad_fn = constraint_grad_fn,
                                                      compute_gradient = TRUE)$gr

                        )
      )


      # do.call(c,
      #
      #         # lapply(1:J,
      #         #   function(j){
      #         #     means <- exp(X%*%B[,j,drop = FALSE] + z);
      #         #   return(-t(X)%*%(Y[,j] - means) +
      #         #   u*(ifelse(j == j_constr,1,0) - con_gr[j])*ek_constr +
      #         #   rho*(B[k_constr,j_constr] -
      #         #   lin_const(B[k_constr,],
      #         #   con_gr,
      #         #   B0[k_constr,]))*
      #         #     (ifelse(j == j_constr,1,0) - con_gr[j])*ek_constr)
      #         #     }
      #         #   )
      # )

      lag_gr_norm <- sqrt(sum(lag_gr^2))
      # print(lag_gr_norm)
      lgns <- c(lgns,lag_gr_norm)

      # print(lag_gr_norm)

    }

    # step <- B - B0
    # step_norm <- max(abs(step))
    #
    # if(step_norm>0.1){
    #   step <- 0.1*step/step_norm
    # }
    # B <- B0 + step

    feas_gap <- B[k_constr,j_constr] - constraint_fn(B[k_constr,])

    if(verbose){
      message("Augmented Lagrangian subproblem solved with u = ", round(u,2), " and rho = ", round(rho,2),
              ". Solution deviates from feasibility by ", signif(feas_gap,3),".")
    }
    u <- u + rho*feas_gap
    # us <- c(us,u)
    rho <- rho*rho_scaling



  }

  return(B)

}

# }
# B_cups %>%
#   group_by(k,j) %>%
#   mutate(value = value - value[iter ==max(iter)]) %>%
#   filter(iter>500) %>%
#   ggplot() +
# 
#   # geom_point(aes(x = iter, y= asinh(value),color = u))
#   geom_line(aes(x = iter,y = asinh(value), group = interaction(k,j)))



#
#
# lagrangian_fit <- function(X,
#                            Y,
#                            B = NULL,
#                            constraint_fn = NULL,
#                            constraint_grad_fn = NULL,
#                            j_constr,
#                            k_constr,
#                            maxit = 250,
#                            tolerance = 1e-1,
#                            info_reg = 0.5,
#                            max_step = 0.5,
#                            c1 = 1e-4,
#                            verbose = TRUE){
#
#   J <- ncol(Y)
#   n <- nrow(Y)
#   p <- ncol(X)
#   X_cup <- X_cup_from_X(X,J)
#   # B[k_constr,j_constr] <- constraint_fn(B[k_constr,-j_constr])
#   B_cup <- B_cup_from_B(B)
#   z_design <- Matrix::kronecker(Matrix::Diagonal(n = n),matrix(1,ncol = 1,nrow = J))
#   z <- update_z_no_wts(Y,X,B)
#   Y_long <- Y_long_from_Y(Y)
#
#
#   g_cup <- function(B_cup,k,p,J){
#     Bk_inds <- k + ((1:J) - 1)*p
#     return(constraint_fn(B_cup[Bk_inds,]))
#   }
#
#   constr_ind <- k_constr + (j_constr - 1)*p
#
#
#   long_means <- exp(X_cup%*%B_cup+ z_design%*%matrix(z,ncol = 1))
#   # u_mult <- B_from_B_cup(Matrix::crossprod(-X_cup,Y_long - long_means),J,p)[k_constr,j_constr]
#   # u <- 0#u_mult#*(B[k_constr,j_constr] - constraint_fn(B[k_constr,]))
#   # rho <- 1e-3*u_mult
#
#   rho <- 100
#   u <- B[k_constr,j_constr] - constraint_fn(B[k_constr,])
#   u <- rho*u
#   us <- u
#   fgs <- B[k_constr,j_constr] - constraint_fn(B[k_constr,])
#
#   # get_constraint_gr <- function(B_cup,k_constr,j_constr){
#   #   Bk_constr_inds <- k_constr + ((1:J) - 1)*p
#   #   long_gr <- 0*B_cup
#   #   long_gr[Bk_constr_inds,] <- -constraint_grad_fn(B_cup[Bk_constr_inds,])
#   #   long_gr[k_constr + (j_constr - 1)*p,] <- 1 + long_gr[k_constr + (j_constr - 1)*p,]
#   #   return(long_gr)
#   # }
#
#
#
#   B_cups <- data.frame(value = B_cup[,1],
#                        j = rep(1:J,each = 2),
#                        k = rep(1:2, J),
#                        iter = 0)
# outer_iter <- 1
#   while(abs(feasibility_gap)>1e-5){
#     constraint_fn_value <- constraint_fn(B[k_constr,])
#     constraint_grad_value <- constraint_grad_fn(B[k_constr,])
#
#
#
#     # long_means <- exp(X_cup%*%B_cup+ z_design%*%matrix(z,ncol = 1))
#     # cg <- get_constraint_gr(B_cup,k_constr,j_constr)
#     # lag_grad <- Matrix::crossprod(-X_cup,Y_long - long_means) +
#     #   u*cg +
#     #   rho*cg*(B[k_constr,j_constr] - constraint_fn(B[k_constr,]))
# plg_norm <- sum(lag_grad^2)
# old_plg_norm <- plg_norm
# iter <- 1
#   while(plg_norm>1e-2){
#
#     for(j in 1:J){
#       lin_starting_vals <- linearized_aug_lag_z()
#     stepsize <- 0.5
#
#
#
#   # # curr_lag <- -sum(Y_long*log(long_means) - long_means) + u*(B_cup[constr_ind,] - g_cup(B_cup,k,p,J))
#   # info <- Matrix::bdiag(lapply(1:J, function(j)t(X)%*%diag(as.numeric(exp(X%*%B[,j,drop = FALSE] + z)))%*%X))
#   # update <- solve(info + 0.01*diag(rep(sqrt(lg_norm),nrow(info))))%*%lag_grad
#   # update[-((j - 1)*p + 1:p)] <- 0
#   # accept <- FALSE
#   # old_plg_norm <- plg_norm
#   while(!accept){
#     prop_B_cup <- B_cup - stepsize*update
#     prop_B <- B_from_B_cup(prop_B_cup,J,p)
#     if(max(abs(prop_B_cup))<=20){
#     prop_long_means <- exp(X_cup%*%prop_B_cup+ z_design%*%matrix(z,ncol = 1))
#     pcg <- get_constraint_gr(prop_B_cup,k_constr,j_constr)
#     prop_lag_grad <- Matrix::crossprod(-X_cup,Y_long - prop_long_means) +
#       u*pcg +
#       rho*pcg*(prop_B[k_constr,j_constr] - constraint_fn(prop_B[k_constr,]))
#     plg_norm <- sum(prop_lag_grad^2)
#     # print(lg_norm/plg_norm)
#     # stepsize <- stepsize/2
#     if(!is.nan(plg_norm)){
#     if(plg_norm<= old_plg_norm){
#       accept <- TRUE
#       # prop_B_cup <- B_cup - stepsize*update
#     }}
#     }
#     stepsize <- 0.5*stepsize
#
#   }
#   message("j = ",j, " norm = ",round(plg_norm,1))
#   if(plg_norm>old_plg_norm){
#     stop()
#   }
#
#   pgs <- c(pgs,plg_norm)
#   B_cup <- prop_B_cup
#
#   B <- B_from_B_cup(B_cup,J,p)
#   for(k in 1:p){
#     B[k,] <- B[k,] - constraint_fn(B[k,])
#   }
#   z <- update_z_no_wts(Y=Y,X = X, B= B)
#
#
#
#
#   }
#   B_cup <- B_cup_from_B(B)
#
#
#
#   # if(plg_norm <=100){
#
#   # u <- u + rho*(B_cup[constr_ind,] - g_cup(B_cup,k,p,J))
#   # }
#
#   B <- B_from_B_cup(B_cup,J,p)
# #
# #   lagrangian_gr <-  constr_lagrangian(Y = Y,
# #                                       X = X,
# #                                       B = B,
# #                                       constraint_fn = constraint_fn,
# #                                       constraint_grad_fn = constraint_grad_fn,
# #                                       k_constr = k_constr,
# #                                       j_constr = j_constr)
# #
# #   print(round(sqrt(sum(lagrangian_gr^2)),3))
#   if(iter %%10 == 0){
#   print(round(plg_norm,2))
#     }
#   iter <- iter + 1
#
#   }
#
# B_cups <- rbind(B_cups,data.frame(value = B_cup[,1],
#                                   j = rep(1:J,each = 2),
#                                   k = rep(1:2, J),
#                                   iter = outer_iter))
# outer_iter <- outer_iter + 1
#
# feasibility_gap <- B_cup[constr_ind,] - g_cup(B_cup,k,p,J)
# u <- u + rho*feasibility_gap
# us <- c(us,u)
# rho <- 1.5*rho
# print(abs(feasibility_gap))
# fgs <- c(fgs,feasibility_gap)
# }
# }
#
# B_cups %>%
#   group_by(k,j) %>%
#   # filter( k == k_constr, j == j_constr) %>%
#   # mutate(value = value - value[iter ==max(iter)]) %>%
#   # filter(iter>3000) %>%
#   ggplot() +
#   geom_line(aes(x = iter, y= value, group = interaction(k,j),color = as.factor(j),
#                 linetype = as.factor(k))) +
#   theme_bw()
# B
# plot(us)
# plot(fgs)
