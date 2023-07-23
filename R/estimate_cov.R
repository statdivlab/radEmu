# estimate_cov <- function(obs,
#                          method = "nlshrink"){
#   if(method == "nlshrink"){
#     filler_obj <- testthat::capture_output(
#     cov_est <- nlshrink::nlshrink_cov(obs,method = "nloptr")
#     )
#   }
#   if(method == "linshrink"){
#     cov_est <- nlshrink::linshrink_cov(obs)
#   }
#   if(method == "plugin"){
#     cov_est <- cov(obs)
# 
#   }
# 
#   if(method == "crossproduct"){
#     n <- nrow(obs)
#     cov_est <- (1/(n +1))*Reduce("+",
#                                   lapply(1:n,
#                                              function(i)
#                                                crossprod(obs[i,,drop = FALSE]))
#                                   )
#   }
#   if(method == "2010RBLW"){
#     cov_est <- CovTools::CovEst.2010RBLW(Umat)$S
#   }
#   if(method == "2010OAS"){
#     cov_est <- CovTools::CovEst.2010OAS(Umat)$S
#   }
# 
#   if(!(method %in% c("nlshrink","plugin","linshrink", "2010RBLW","2010OAS","crossproduct"))){
#     stop("Method must be one of nlshrink,linshrink,2010RBLW, 2010OAS, crossproduct, or plugin.")
#   }
# 
#   return(cov_est)
# }
