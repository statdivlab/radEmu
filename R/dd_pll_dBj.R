#
#
# dd_pll_dBj <- function(X,
#                        Y,
#                        B,
#                        j){
#   z <- update_z(Y,X,B)
#
#   ddBj <- numeric(nrow(B))
#   n <- nrow(Y)
#   J <- ncol(Y)
#
#   for(i in 1:n){
#     pij <- exp(X[i,,drop= FALSE]%*%B[,j,drop = FALSE])/
#       sum(exp(X[i,,drop= FALSE]%*%B))
#     for(jprime in 1:J){
#       indicator <- as.numeric(jprime == j)
#       ddBj <- ddBj +
#         as.numeric((Y[i,jprime] - exp(X[i,,drop = FALSE]%*%B[,jprime,drop = FALSE]))*(indicator - pij))*t(X[i,,drop = FALSE])
#     }
#   }
#   return(ddBj)
# }
