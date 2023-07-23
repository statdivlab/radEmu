

info_ljaug <- function(Y,
                           X,
                           B,
                           z,
                           constraint_fn,
                           constraint_grad_fn,
                           j,
                           j_constr,
                           k_constr
){

  mu_j <- exp(X%*%B[,j,drop = FALSE] + z)

  W_j <- diag(as.numeric(mu_j))

  info_j <- t(X)%*%W_j%*%X

  mu_j_constr <- exp(X%*%B[,j_constr,drop = FALSE] + z)
  W_j_constr <- diag(as.numeric(mu_j_constr))

  info_aug_scalar <- as.numeric(
    t(X[,k_constr,drop = FALSE])%*%W_j_constr%*%X[,k_constr,drop = FALSE]
  )

  e_k_constr <- matrix(as.numeric((1:(ncol(X))) == k_constr),ncol = 1)

  dg <- constraint_grad_fn(B[k_constr,-j_constr])
  dg_index <- ifelse(j>j_constr,j - 1,j)
  dg <- dg[dg_index]
  info_aug <- (e_k_constr%*%t(e_k_constr))*info_aug_scalar*dg^2

  return(info_j + info_aug)

}
