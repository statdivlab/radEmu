# emuFit <-   function(X,
#                        Y,
#                        B = NULL,
#                        constraint_fn = huber_center,
#                        maxit = 500,
#                        tolerance = 1e-5,
#                        penalize = TRUE){
#
#   J <- ncol(Y)
#   p <- ncol(X)
#
#   if(penalize){
#     full_fit <- emuFit_micro_penalized(X,
#                                        Y,
#                                        B,
#                                        constraint_fn,
#                                        maxit = maxit,
#                                        tolerance = tolerance)
#
#     for(k in 2:p){
#
#       for(j in 1:J){
#         B_start <- full_fit$B
#         B_start[k,j] <- huber_center(B_start[k,])
#       micro_score_test(Y = full_fit$Y_augmented,
#                        X = X,
#                        B =  B_start,
#                        constraint_fn = huber_center,
#                        constraint_type = "huber",
#                        null_k = k,
#                        null_j = j,
#                        maxit = 500,
#                        huber_param = 1,
#                        tolerance = 1e-2)
#         }
#     }
#
#   }
#
#
#
# }
