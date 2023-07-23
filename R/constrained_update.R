#
#
# constrained_update <-
#   function(B,
#            z,
#            j,
#            k_constr,
#            j_constr,
#            constraint_fn){
#     if(j != j_constr){
#
#       mu_j <- exp(z + X%*%B[,j,drop = FALSE])
#       lj <- sum(Y[,j]*log(mu_j) - mu_j)
#
#       mu_j_constr <- exp(z + X%*%B[,j_constr,drop = FALSE])
#       lj_constr <-  sum(Y[,j_constr]*log(mu_j_constr) - mu_j_constr)
#
#       l_aug <- lj + lj_constr
#
#       dlj_dBj <- t(X)%*%(Y[,j] - mu_j) #gradient contribution of Y[,j]
#       dlj_constr_dBj_constr <- t(X)%*%(Y[,j_constr] - mu_j_constr)
#
#       cf <- function(Bkconstrj){
#         Bkconstr <- B[k_constr,]
#         Bkconstr[j] <- Bkconstrj
#         constraint_fn(Bkconstr[-j_constr])
#       }
#
#       dh <- numDeriv::grad(cf,B[k_constr,j])
#
#       for_dh_analytic <- as.numeric(B[k_constr,] - huber_center(B[k_constr,-j_constr])<= huber_param)
#       dh_analytic <- 2*for_dh_analytic[j]/sum(for_dh_analytic[-j_constr])
#
#       dh_analytic/dh
#
#       dlj_constr_dBjk_constr <- dlj_constr_dBj_constr[k_constr,]*dh
#
#       dlaug_dBj <- dlj_dBj + dlj_constr_dBjk_constr*as.numeric(1:p == k_constr)
#
#       ljaug <- function(Bj){
#         Bj <- matrix(Bj,ncol = 1) #input value of Bj
#         lj <- sum(Y[,j]*(X%*%Bj + z) - exp(X%*%Bj +z)) #mean for Yj
#         Bjconstr <- B[,j_constr,drop = FALSE]
#         Bkconstr <- B[k_constr,]
#         Bkconstr[j] <- Bj[k_constr] # update with k_constr val of Bj
#         Bjconstr <- B[,j_constr,drop = FALSE]
#         Bjconstr[k_constr,] <- constraint_fn(Bkconstr[-j_constr])
#         aug <- sum(Y[,j_constr]*(X%*%Bjconstr + z) - exp(X%*%Bjconstr +z))
#         return(lj + aug)
#       }
#
#       numerical_deriv <- numDeriv::grad(ljaug,B[,j])
#       numerical_hess <- numDeriv::hessian(ljaug,B[,j])
#
#       info_j <- t(X)%*%diag(as.numeric(mu_j))%*%X
#       info_j[k_constr,k_constr] <-
#         info_j[k_constr,k_constr] +  (t(X)%*%diag(as.numeric(mu_j_constr))%*%X)[k_constr,k_constr]*dh^2
#
#
#       update <- 0.5*qr.solve(info_j + 0.001*diag(as.numeric(dlaug_dBj^2)), dlaug_dBj,tol = 1e-20)
#
#
#       l_aug_prop <- l_aug - 1
#
#       while(l_aug_prop < l_aug){
#       prop_B <- B
#       prop_B[,j] <- prop_B[,j] + update
#       prop_B[k_constr,j_constr] <- constraint_fn(prop_B[k_constr,-j_constr])
#
#       prop_mu_j <- exp(z + X%*%prop_B[,j,drop = FALSE])
#       prop_lj <- sum(Y[,j]*log(prop_mu_j) - prop_mu_j)
#
#       prop_mu_j_constr <- exp(z + X%*%prop_B[,j_constr,drop = FALSE])
#       prop_lj_constr <-  sum(Y[,j_constr]*log(prop_mu_j_constr) - prop_mu_j_constr)
#
#       l_aug_prop<- prop_lj + prop_lj_constr
#
#       update <- update/2
#
#       }
#       return(prop_B[,c(j,j_constr)])
#     }
#
#     if(j == j_constr){
#       mu_j <- exp(z + X%*%B[,j,drop = FALSE])
#       lj <- sum(Y[,j]*log(mu_j) - mu_j)
#
#       dlj_dBj <- t(X)%*%(Y[,j] - mu_j)
#       info_j <- t(X)%*%diag(as.numeric(mu_j))%*%X
#       dlj_dBj <- dlj_dBj[-k_constr,drop = FALSE]
#       info_j <- info_j[-k_constr,-k_constr,drop = FALSE]
#
#       update <- B[,j,drop = FALSE]*0
#       update[-k_constr,] <- 0.5*qr.solve(info_j,dlj_dBj)
#
#       lj_prop <- lj - 1
#
#       while(lj_prop < lj){
#         prop_B <- B
#         prop_B[,j] <- prop_B[,j] + update
#
#         prop_mu_j <- exp(z + X%*%prop_B[,j,drop = FALSE])
#         lj_prop <- sum(Y[,j]*log(prop_mu_j) - prop_mu_j)
#
#         update <- update/2
#
#       }
#       return(prop_B[,j])
#
#     }
#   }
