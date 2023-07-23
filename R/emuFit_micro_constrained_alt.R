#
# emuFit_micro_constrained_alt <-
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
#            # step_decay = TRUE,
#            # abs_max_step = 0.5,
#            c1 = 1e-4,
#            step_ratio = 0.5,
#            verbose = TRUE
#   ){
#     n <- nrow(Y)
#     J <- ncol(Y)
#     p <- ncol(X)
#     orig_step <- max_step
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
#
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
#     lls <- lgs <-  numeric(0)
#
#     # B_list <- vector(0,mode = "list")
#     B_cup <- B_cup_from_B(B)
#
#     B_cup_df <- data.frame(matrix(ncol= nrow(B_cup) , nrow = maxit +1))
#     B_cup_df$iter <- 0:maxit
#     B_cup_df[1,1:(p*J)] <- as.numeric(as.matrix(B_cup))
#
#
#
#
#     while(!converged){
#
#       # long_means <- exp(X_cup%*%B_cup + z_design%*%matrix(z,ncol = 1))
#
#       means <- exp(X%*%B + matrix(z,ncol = 1)%*%matrix(1,nrow = 1, ncol = J))
#
#       neg_ll_B_grad <- -do.call(cbind,lapply(1:J,function(j) t(X)%*%(
#         Y[,j,drop = FALSE] - means[,j,drop = FALSE]
#         )
#         ))
#       neg_ll <-  sum(Y*log(means) - prop_means)
#
#
#
#       stepsize <- max_step
#       accept <- FALSE
#       while(!accept){
#
#         prop_B <- B + stepsize*neg_ll_B_grad
#         prop_B[k_constr,j_constr] <- constraint_fn(prop_B[k_constr,-j_constr])
#
#         update <- prop_B - B
#
#         prop_means <- exp(X%*%prop_B + matrix(z,ncol = 1)%*%matrix(1,nrow = 1, ncol = J))
#
#        neg_prop_ll <- -sum(Y*log(prop_means) - prop_means)
#
#        if(is.nan(neg_prop_ll)){
#          accept <- FALSE
#        } else{
#        if(neg_prop_ll <= neg_ll + c1*sum(neg_ll_B_grad*update)){
#          accept <- TRUE
#        }
#        }
#
#        stepsize <- step_ratio*stepsize
#       }
#
#       lls <- c(lls,-neg_ll)
#
#       B <- prop_B
#       for(k in 1:p){
#         B[k,] <- B[k,] - constraint_fn(B[k,])
#       }
#
#       z <- update_z_no_wts(Y,X,B)
#
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
#       #
#       # deriv_norm <- sqrt(sum(lagrangian_gr^2))
#       # nlg <- length(lgs)
#       # n_to_plot <- min(nlg,100)
#       # if(lgs[nlg] - lgs[nlg - 1]>0){
#       #   warning("norm of gradient increased!")
#       # }
#       #
#       # plot(1:n_to_plot + nlg - n_to_plot,lgs[(nlg - n_to_plot + 1):nlg],type = "l")
#
#
#       if(verbose){
#         message("Norm of lagrangian derivative ",signif(deriv_norm,4) )
#       }
#       print(iter)
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
#
