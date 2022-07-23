# test_that("Block inversion of information matrix returns correct result", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each =4))
#   z <- rnorm(8) +8
#   b0 <- rnorm(10)
#   b1 <- 0*b0
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 10, nrow = 8)
#   for(i in 1:8){
#     for(j in 1:10){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rpois(1, lambda = temp_mean)
#     }
#   }
#
#   p <- ncol(X)
#   J <- ncol(Y)
#   n <- nrow(X)
#
#
#   B <- b
#
#
#
#   weights <- rep(1,n*J)
#   rect_weights <- 1.0 + 0.0*Y
#
#     X_tilde <- X_to_X_tilde(X,J)
#     # Y_tilde <- Y_to_Y_tilde(Y)
#     S <- Matrix::sparseMatrix(i = 1:(n*J),
#                               j = rep(1:n,each = J),
#                               x = rep(1, n*J))
#     D_tilde <- cbind(X_tilde,S)
#
#     z <- apply(Y*rect_weights,1,function(x) log(sum(x))) -
#       apply(exp(X%*%B)*rect_weights,1,function(x) log(sum(x)))
#
#     rownames(D_tilde) <- 1:(n*J)
#     B_tilde <- B_to_B_tilde(B)
#     theta <- rbind(B_tilde,Matrix::Matrix(z,ncol = 1))
#     X_tilde_repar <- X_tilde
#     X_tilde_repar <- X_tilde_repar[,-((1:p)*J)]
#     D_tilde_repar <- cbind(X_tilde_repar,S)
#
#     B_tilde <- B_to_B_tilde(B)
#     theta <- rbind(B_tilde,Matrix::Matrix(z,ncol = 1))
#     W <- Matrix::Diagonal(x = weights*as.numeric(exp((D_tilde%*%theta))))
#
#     B_tilde <- B_to_B_tilde(B)
#     theta <- rbind(B_tilde,Matrix::Matrix(z,ncol = 1))
#     info <- Matrix::crossprod(D_tilde_repar,W)%*%D_tilde_repar
#     info_inv_direct <-
#         Matrix::solve(info)
#
#     info_inv_block <-
#       get_info_inv(X_tilde = X_tilde_repar,
#                    S = S,
#                    n = n,
#                    J = J,
#                    glm_weights = weights*as.numeric(exp((D_tilde%*%theta))))
#
#     expect_equal(as.matrix(info_inv_block),as.matrix(info_inv_direct))
#
#
#
#
#
#
# })
