

micro_fisher <- function(X,Yj,Bj,z,
                         stepsize = 1,
                         c1 = 0.1){
  means <- exp(X%*%Bj + z)
  info <- t(X)%*%diag(as.numeric(means))%*%X
  lj_grad <- colSums(diag(as.numeric(Yj - means))%*%X)

  update <- try(stop(),silent = TRUE)


  if(nrow(info) >1){
    info_avg_diag <- diag(rep(sqrt(mean(diag(info)^2)),nrow(info)))
  } else{
    info_avg_diag <- abs(info)
  }

  regularization <- 0
  while(inherits(update,"try-error")){
    update <- try(qr.solve(info + regularization*info_avg_diag,lj_grad),silent = TRUE)
    regularization <- ifelse(regularization ==0, 0.01,10*regularization)
  }


  obj <- -sum(Yj*log(means) - means)
  obj_grad <- -lj_grad

  suff_decrease_term <- c1*sum(obj_grad*update)

  suff_decrease <- FALSE
  while(!(suff_decrease)){

    prop_Bj <- Bj + stepsize*update
    prop_log_mu <- X%*%prop_Bj + z
    prop_obj <- -sum(Yj*prop_log_mu - exp(prop_log_mu))

    suff_decrease <- prop_obj <= obj + suff_decrease_term*stepsize

    stepsize <- stepsize*0.5

  }

  stepsize <- stepsize/0.5

  return(stepsize*update)
}


#
#
# micro_fisher <- function(X,Yj,Bj,z,
#                          stepsize = 0.1,
#                          c1 = 1e-4,
#                          step_ratio = 0.5){
#   means <- exp(X%*%Bj + z)
#   info <- t(X)%*%diag(as.numeric(means))%*%X
#   lj_grad <- colSums(diag(as.numeric(Yj - means))%*%X)
#
#   update <- try(stop(),silent = TRUE)
#
#
#   if(nrow(info) >1){
#   info_avg_diag <- diag(rep(sqrt(mean(diag(info)^2)),nrow(info)))
#   } else{
#     info_avg_diag <- abs(info)
#   }
#
#   regularization <- 0
#   while(inherits(update,"try-error")){
#     update <- try(qr.solve(info + regularization*info_avg_diag,lj_grad),silent = TRUE)
#     regularization <- ifelse(regularization ==0, 0.01,10*regularization)
#   }
#
#   # if(max(abs(update))> abs_max_step){
#   #   update <- abs_max_step*update/max(abs(update))
#   # }
#
#
#
#   obj <- -sum(Yj*log(means) - means)
#   obj_grad <- -lj_grad
#
#   suff_decrease_term <- c1*sum(obj_grad*update)
#
#   curv_term <- c2*sum(obj_grad*update)
#
#   suff_decrease <- FALSE
#   curv_cond <- FALSE
#
#   ##fiddling
#   # eps <- exp(seq(0,-10,-0.1))
#   # obj_fn <- function(x){
#   #   temp_Bj <- Bj + x*update
#   #   temp_means <- exp(X%*%temp_Bj + z)
#   #   return(-sum(Yj*log(temp_means) - temp_means))
#   # }
#   #
#   # gr_fn <- function(x){
#   #   temp_Bj <- Bj + x*update
#   #   temp_means <- exp(X%*%temp_Bj + z)
#   #   return(-t(X)%*%(Yj - temp_means))}
#   # objs <- sapply(eps, function(x) obj_fn(x))
#   # grs <- sapply(eps, function(x) gr_fn(x))
#
#
#   while( (!suff_decrease)|(!curv_cond)){
#
#     prop_Bj <- Bj + stepsize*update
#     prop_log_mu <- X%*%prop_Bj + z
#     prop_obj <- -sum(Yj*prop_log_mu - exp(prop_log_mu))
#
#     suff_decrease <- prop_obj <= obj + suff_decrease_term*stepsize
#
#     # if(suff_decrease){
#     #   prop_means <-  exp(X%*%prop_Bj + z)
#     #   prop_grad <- -1*colSums(diag(as.numeric(Yj - prop_means))%*%X)
#     #   curv_cond <- sum(prop_grad*update) <= curv_term
#     # }
#     curv_cond <- TRUE
#
#     stepsize <- stepsize*step_ratio
#
#   }
#
#   stepsize <- stepsize/0.95
#
#   return(stepsize*update)
# }
#
#
#
#
# #
# # micro_fisher <- function(X,Yj,Bj,z,regularization = 0.5){
# #   means <- exp(X%*%Bj + z)
# #   info <- t(X)%*%diag(as.numeric(means))%*%X
# #   lj_grad <- colSums(diag(as.numeric(Yj - means))%*%X)
# #   return(qr.solve(info + regularization*sqrt(sum(lj_grad^2))*diag(nrow(info)),lj_grad,
# #                   tol = 1e-20))
# # }
