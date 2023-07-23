
micro_fisher_constr <- function(X,Y,B,z,j,j_constr,
                                k_constr,
                                constraint_fn,
                                constraint_grad_fn,
                                c1 = 0.1,
                                stepsize= 0.5,
                                abs_max_step = 0.5){

  if(j != j_constr){
    aug_grad <- dljaug_dBj(Y = Y,
                           X = X,
                           B = B,
                           z = z,
                           constraint_fn = constraint_fn,
                           constraint_grad_fn = constraint_grad_fn,
                           j = j,
                           j_constr = j_constr,
                           k_constr = k_constr)
    aug_info <- info_ljaug(Y = Y,
                           X = X,
                           B = B,
                           z = z,
                           constraint_fn = constraint_fn,
                           constraint_grad_fn = constraint_grad_fn,
                           j = j,
                           j_constr = j_constr,
                           k_constr = k_constr)

    update <- try(stop(),silent = TRUE)

    aug_diag <- sqrt(mean(diag(aug_info)^2))
    aug_diag <- diag(rep(1,nrow(aug_info)))*aug_diag
    regularization <- 0
    while(inherits(update,"try-error")){

      update <- try(qr.solve(aug_info + regularization*aug_diag,
                             aug_grad),silent = TRUE)

      regularization <- ifelse(regularization ==0, 0.1,10*regularization)
    }

    if(max(abs(update))>abs_max_step){
      update <- abs_max_step*update/max(abs(update))
    }

    ### line search

    neg_aug_ll_j <- function(Bj,temp_B){
      temp_B[,j] <- Bj
      temp_B[k_constr,j_constr] <- constraint_fn(temp_B[k_constr,-j_constr])

      log_mu <- X%*%temp_B[,c(j,j_constr)] + cbind(z,z)

      return(-sum(Y[,c(j,j_constr)]*log_mu - exp(log_mu)))
    }

    obj <-  neg_aug_ll_j(B[,j],B)

    obj_grad <- -aug_grad

    suff_decrease_term <- c1*sum(obj_grad*update)
    curv_term <- abs(sum(obj_grad*update))

    suff_decrease <- FALSE
    curv_cond <- FALSE
    while(!(suff_decrease)){


      prop_Bj <- B[,j] + stepsize*update
      prop_obj <- neg_aug_ll_j(prop_Bj,B)

      suff_decrease <- prop_obj <= obj + suff_decrease_term*stepsize

      # if(suff_decrease){
      #   prop_B <- B
      #   prop_B[,j] <- prop_Bj
      #   prop_B[k_constr,j_constr] <- constraint_fn(B[k_constr,-j_constr])
      #   prop_grad <- -dljaug_dBj(Y = Y,
      #                           X = X,
      #                           B = prop_B,
      #                           z = z,
      #                           constraint_fn = constraint_fn,
      #                           constraint_grad_fn = constraint_grad_fn,
      #                           j = j,
      #                           j_constr = j_constr,
      #                           k_constr = k_constr)
      #
      #   curv_cond <- abs(sum(prop_grad*update)) <= curv_term
      # } else{
      #   curv_cond <- FALSE
      # }

      stepsize <- stepsize*0.5
      # print(stepsize)

    }
    stepsize <- stepsize/0.5
    return(update*stepsize)






  } else{
    partial_update <- micro_fisher(X = X[,-k_constr,drop = FALSE],
                                   Yj = Y[,j, drop = FALSE],
                                   Bj = B[-k_constr,j,drop = FALSE],
                                   z = z,
                                   c1 = c1,
                                   stepsize = stepsize)
    if(max(abs(partial_update))>abs_max_step){
      partial_update <- abs_max_step*partial_update/max(abs(partial_update))
    }
    update <- 0*B[,j]
    update[-k_constr] <- partial_update
    return(update)
  }



}
# micro_fisher_constr <- function(X,Y,B,z,j,j_constr,
#                                 k_constr,
#                                 constraint_fn,
#                                 constraint_grad_fn,
#                                 c1 = 1e-4,
#                                 step_ratio = 0.5,
#                                 # abs_max_step = 1,
#                                 stepsize= 0.1){
#
#   if(j != j_constr){
#     aug_grad <- dljaug_dBj(Y = Y,
#                            X = X,
#                            B = B,
#                            z = z,
#                            constraint_fn = constraint_fn,
#                            constraint_grad_fn = constraint_grad_fn,
#                            j = j,
#                            j_constr = j_constr,
#                            k_constr = k_constr)
#     aug_info <- info_ljaug(Y = Y,
#                            X = X,
#                            B = B,
#                            z = z,
#                            constraint_fn = constraint_fn,
#                            constraint_grad_fn = constraint_grad_fn,
#                            j = j,
#                            j_constr = j_constr,
#                            k_constr = k_constr)
#
#     update <- try(stop(),silent = TRUE)
#
#     aug_diag <- sqrt(mean(diag(aug_info)^2))
#     aug_diag <- diag(rep(1,nrow(aug_info)))*aug_diag
#     regularization <- 0
#     while(inherits(update,"try-error")){
#
#     update <- try(qr.solve(aug_info + regularization*aug_diag,
#                            aug_grad),silent = TRUE)
#
#     regularization <- ifelse(regularization ==0, 0.1,10*regularization)
#     }
#
#     # if(max(abs(update))> abs_max_step){
#     #   update <- abs_max_step*update/max(abs(update))
#     # }
#
#     ### line search
#
#     neg_aug_ll_j <- function(Bj,temp_B){
#       temp_B[,j] <- Bj
#       temp_B[k_constr,j_constr] <- constraint_fn(temp_B[k_constr,-j_constr])
#
#       log_mu <- X%*%temp_B[,c(j,j_constr)] + cbind(z,z)
#
#       return(-sum(Y[,c(j,j_constr)]*log_mu - exp(log_mu)))
#     }
#
#     ## fiddling
#     # eps <- exp((seq(0,-10,-0.01)))
#     # objs <- sapply(eps,function(e) neg_aug_ll_j(B[,j] + e*update,B))
#     #
#     # aug_grad_fn <- function(Bj,temp_B){
#     #   temp_B[,j] <- Bj
#     #   temp_B[k_constr,j_constr] <- constraint_fn(temp_B[k_constr,-j_constr])
#     #   return(abs(sum(update*dljaug_dBj(Y = Y,
#     #                        X = X,
#     #                        B = temp_B,
#     #                        z = z,
#     #                        constraint_fn = constraint_fn,
#     #                        constraint_grad_fn = constraint_grad_fn,
#     #                        j = j,
#     #                        j_constr = j_constr,
#     #                        k_constr = k_constr))))
#     # }
#     #
#     # grnorms <-  sapply(eps,function(e)     aug_grad_fn(B[,j] + e*update,B))
#     #
#     # plot(eps,grnorms)
#     # plot(eps,objs)
#
#     obj <-  neg_aug_ll_j(B[,j],B)
#
#     obj_grad <- -aug_grad
#
#     suff_decrease_term <- c1*sum(obj_grad*update)
#     curv_term <- c2*sum(obj_grad*update)
#
#     suff_decrease <- FALSE
#     curv_cond <- FALSE
#     while( (!suff_decrease)|(!curv_cond)){
#
#
#     prop_Bj <- B[,j] + stepsize*update
#     prop_obj <- neg_aug_ll_j(prop_Bj,B)
#
#     suff_decrease <- prop_obj <= obj + suff_decrease_term*stepsize
#
#     # if(suff_decrease){
#     #   prop_B <- B
#     #   prop_B[,j] <- prop_Bj
#     #   prop_B[k_constr,j_constr] <- constraint_fn(B[k_constr,-j_constr])
#     #   prop_grad <- -dljaug_dBj(Y = Y,
#     #                           X = X,
#     #                           B = prop_B,
#     #                           z = z,
#     #                           constraint_fn = constraint_fn,
#     #                           constraint_grad_fn = constraint_grad_fn,
#     #                           j = j,
#     #                           j_constr = j_constr,
#     #                           k_constr = k_constr)
#     #
#     #   curv_cond <- sum(prop_grad*update) <= curv_term
#     # } else{
#     #   curv_cond <- FALSE
#     # }
#     curv_cond <- TRUE
#
#     stepsize <- stepsize*step_ratio
#     # print(stepsize)
#
#     }
#     return(update*stepsize)
#
#
#
#
#
#
#   } else{
#     partial_update <- micro_fisher(X = X[,-k_constr,drop = FALSE],
#                                    Yj = Y[,j, drop = FALSE],
#                                    Bj = B[-k_constr,j,drop = FALSE],
#                                    z = z,
#                                    stepsize = max_step,
#                                    step_ratio = step_ratio)
#     update <- 0*B[,j]
#     update[-k_constr] <- partial_update
#     return(update)
#   }
#
#
#
# }
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# # micro_fisher_constr <- function(X,Y,B,z,j,jstar,kstar,dg_dBj){
# #
# #   gr = constr_grad(X,Y,B,z,j,jstar,kstar,dg_dBj)
# #   info = constr_info(X,Y,B,z,j,jstar,kstar,dg_dBj)
# #
# #   dampening <- 0.5
# #   update <- suppressWarnings(try(stop(),silent = TRUE))
# #   if(j != jstar){
# #     while(inherits(update,"try-error")){
# #     update <- try(as.numeric(qr.solve(info + dampening*sqrt(sum(gr^2))*diag(nrow(info)),t(gr))),silent = TRUE)
# #     dampening <- dampening*2
# #     }
# #     return(update)
# #   } else{
# #     gr <- gr[,-kstar,drop = FALSE]
# #     info <- info[-kstar,-kstar,drop = FALSE]
# #     while(inherits(update,"try-error")){
# #       update <- try(as.numeric(qr.solve(info + dampening*sqrt(sum(gr^2))*diag(nrow(info)),t(gr))),silent = TRUE)
# #       dampening <- dampening*2
# #     }
# #     new_Bj <- rep(0,nrow(B))
# #     new_Bj[-kstar] <- update
# #
# #     return(new_Bj)
# #   }
# #
# #
# # }
