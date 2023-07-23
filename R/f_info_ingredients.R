f_info_ingredients <- function(Y,
                   B_cup,
                   B,
                   X,
                   X_cup){
  n <- nrow(Y)
  J <- ncol(Y)
  p <- ncol(X)

  f_info <- matrix(0,ncol = p*J,nrow = p*J)

  z <- update_z_no_wts(Y,X,B)

  M_is <- vector(n,mode = "list")

  for(i in 1:n){
    which_rows <- 1:J + (i - 1)*J
    mu_i <- exp(X_cup[which_rows,]%*% B_cup + z[i])
    p_i <- mu_i/sum(mu_i)

    M_is[[i]] <- diag(rep(1,J)) - matrix(1,nrow=J,ncol=1)%*%Matrix::t(p_i)
    M_is[[i]] <- Matrix::t(M_is[[i]])%*%diag(as.numeric(mu_i))%*%M_is[[i]]
  }

  W <- Matrix::bdiag(M_is)

  return(W)

}
