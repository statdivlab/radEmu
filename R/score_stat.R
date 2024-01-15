# score_stat <- function(X_cup,
#                        X,
#                        Y,
#                        B,
#                        constraint_grad_fn,
#                        indexes_to_remove,
#                        n,
#                        p,
#                        J){
#   for(i in 1:n){
#     X_cup_i <- X_cup[(i - 1)*J + 1:J,]
#     scores[[i]] <- as.matrix(dpll_dB_cup(X[i,,drop = FALSE],Y[i,,drop = FALSE],B))
#   }
#
#   H <- matrix(0,nrow = p, ncol = J)
#   for(j in 1:J){
#     # H[k_constr,j] <- get_dg_dBj(B = B,
#     #                           j = j,
#     #                           k_star = k_constr,
#     #                           huber_param = huber_param,
#     #                           constraint_type = constraint_type,
#     #                           for_testing = TRUE)[k_constr]
#
#     H[k_constr,j] <- constraint_grad_fn(B[k_constr,])[j]
#   }
#
#
#   H_cup <- B_cup_from_B(H)
#
#   H_cup <- H_cup[-indexes_to_remove,,drop = FALSE]
#
#   B_cup <- B_cup_from_B(B)
#   I <- f_info(Y,B_cup,B,X,X_cup)
#   I <- I[-indexes_to_remove,-indexes_to_remove]
#   # I_inv <- solve(I)
#   # H_I_inv <- crossprod(H_cup,I_inv)
#   H_I_inv <- try(t(qr.solve(I,H_cup,tol= 1e-20)))
#   if(inherits(H_I_inv,"try-error")){
#     H_I_inv <- t(MASS::ginv(as.matrix(I))%*%H_cup)
#     warning("Information matrix numerical singular; check convergence of fit under null.")
#   }
#
#   Dy <- Reduce("+",lapply(scores,function(x) tcrossprod(x)))
#   Dy <- Dy[-indexes_to_remove,-indexes_to_remove]
#   score <- Reduce("+",scores)
#   score <- as.matrix(score)
#   score <- score[-indexes_to_remove,,drop =FALSE]
#
#
#   outside <- H_I_inv%*%score
#
#   inside <- H_I_inv%*%Dy%*%t(H_I_inv)
#
#   score_stat <- (n/(n - 1))*
#     as.numeric(outside)^2/as.numeric(inside)
#   return(score_stat)
# }
