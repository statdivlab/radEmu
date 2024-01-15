# micro_fisher_null <- function(X,
#                               Y,
#                               B,
#                               j,
#                               u,
#                               rho,
#                               k_constr,
#                               constraint_fn,
#                               constraint_grad_fn,
#                               constraint_hess_fn,
#                               z,
#                          stepsize = 0.5,
#                          max_step = 1,
#                          c1 = 0.1){
#
#   means <- exp(X%*%B[,j,drop = FALSE] + z)
#
#   cf <- constraint_fn(B[k_constr,])
#   cg <- constraint_grad_fn(B[k_constr,])
#   ch <- constraint_hess_fn(B[k_constr,],j,j)
#
#   # num_cg <- numDeriv::grad(constraint_fn,B[k_constr,])
#   #
#   # num_ch <- numDeriv::grad(function(x){
#   #   tb <- B[k_constr,];
#   #   tb[j] <- x;
#   #   constraint_fn(tb)},
#   #   B[k_constr,j])
#
#   obj_grad <- -crossprod(X,Y[,j] - means)
#   obj_grad[k_constr,] <- obj_grad[k_constr,] +
#     (u + rho*cf)*cg[j]
# #
# #   obj_numgrad <- numDeriv::grad(function(x){
# #     log_means <- X%*%matrix(x,ncol = 1) + z;
# #     tbk_constr <- B[k_constr,];
# #     tbk_constr[j] <- x[k_constr];
# #     tcf <- constraint_fn(tbk_constr);
# #     return(-sum(Y[,j]*log_means - exp(log_means)) + u*tcf +
# #       (rho/2)*tcf^2)},
# #     B[,j])
#
#
#   obj_hess <- crossprod(X,diag(as.numeric(means)))%*%X
#   obj_hess[k_constr,k_constr] <-
#     obj_hess[k_constr,k_constr] + rho*(cg[j])^2 +
#     ch*(u + rho*cf)
#
#   # obj_numhess <- numDeriv::hessian(function(x){
#   #   log_means <- X%*%matrix(x,ncol = 1) + z;
#   #   tbk_constr <- B[k_constr,];
#   #   tbk_constr[j] <- x[k_constr];
#   #   tcf <- constraint_fn(tbk_constr);
#   #   -sum(Y[,j]*log_means - exp(log_means)) + u*tcf +
#   #      (rho/2)*tcf^2
#   #   },
#   #   B[,j],
#   #   method = "complex")
#   #
#   # # message("Hessian deviates from numerical approx by ",signif(sum(abs(obj_hess - obj_numhess)),3),".")
#   # abs_hess <- sum(abs(obj_hess - obj_numhess))
#   # prop_hess <- max(abs(obj_hess/obj_numhess - 1))
#   #
#   # if(abs_hess >1 & prop_hess >0.5){
#   #   stop("Oops looks like a problem with the hessian")
#   # }
# #
# #   abs_grad <- sum(abs(obj_grad - obj_numgrad))
# #   prop_grad <- max(abs(obj_grad/obj_numgrad - 1))
# #
# #   if(abs_grad >1 & prop_grad >0.5){
# #     stop("Oops looks like a problem with the gradient")
# #   }
#
#
#
#
#
#   update <- try(stop(),silent = TRUE)
#
#
#   if(nrow(obj_hess) >1){
#     hess_avg_diag <- diag(rep(sqrt(mean(diag(obj_hess )^2)),nrow(obj_hess )))
#   } else{
#     hess_avg_diag <- abs(obj_hess )
#   }
#
#   regularization <- 0
#   while(inherits(update,"try-error")){
#     update <- try(-qr.solve(obj_hess + regularization*hess_avg_diag,obj_grad),silent = TRUE)
#     # update <- try(-qr.solve(obj_hess + regularization*diag(abs(diag(obj_hess))),obj_grad),silent = TRUE)
#     regularization <- ifelse(regularization ==0, 0.01,10*regularization)
#   }
#
#   obj <- -sum(Y[,j]*log(means) - means) + u*cf +
#     (rho/2)*cf^2
#
#
#   suff_decrease_term <- c1*sum(obj_grad*update)
#
#   while(suff_decrease_term >= 0){
#
#     update <- try(-qr.solve(obj_hess + regularization*hess_avg_diag,obj_grad),silent = TRUE)
#     # update <- try(-qr.solve(obj_hess + regularization*diag(abs(diag(obj_hess))),obj_grad),silent = TRUE)
#     regularization <- ifelse(regularization ==0, 0.01,10*regularization)
#     suff_decrease_term <- c1*sum(obj_grad*update)
#   }
#
#   if(max(abs(update))>max_step){
#     update <- max_step*update/max(abs(update))
#     suff_decrease_term <- c1*sum(obj_grad*update)
#   }
#
#   suff_decrease <- FALSE
#   while(!(suff_decrease)){
#
#     prop_Bj <- B[,j] + stepsize*update
#     prop_Bk_constr <- B[k_constr,]
#     prop_Bk_constr[j] <- prop_Bj[k_constr]
#     prop_log_mu <- X%*%prop_Bj + z
#     prop_cf <- constraint_fn(prop_Bk_constr)
#     prop_obj <- -sum(Y[,j]*prop_log_mu - exp(prop_log_mu)) +
#       u*prop_cf + (rho/2)*prop_cf^2
#
#     suff_decrease <- prop_obj <= obj + suff_decrease_term*stepsize
#
#     stepsize <- stepsize*0.8
#
#   }
#
#   stepsize <- stepsize/0.5
#
#   return(list(update = stepsize*update,updated_gap = prop_cf))
# }
#
#
