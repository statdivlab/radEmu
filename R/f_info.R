f_info <- function(Y,
                   B_cup,
                   B,
                   X,
                   X_cup){
  n <- nrow(Y)
  J <- ncol(Y)
  p <- ncol(X)
  X_cup <- X_cup_from_X(X,J)

  f_info <- matrix(0,ncol = p*J,nrow = p*J)

  z <- update_z_no_wts(Matrix::Matrix(Y),Matrix::Matrix(X),B)

  for(i in 1:n){
    which_rows <- 1:J + (i - 1)*J
    mu_i <- exp(X_cup[which_rows,]%*% B_cup + z[i])
    p_i <- mu_i/sum(mu_i)

    M_iX_i <- (diag(rep(1,J)) - matrix(1,nrow=J,ncol=1)%*%Matrix::t(p_i))%*%X_cup[which_rows,]

    f_info <- f_info + Matrix::t(M_iX_i)%*%Matrix::Diagonal(x= as.numeric(mu_i))%*%M_iX_i
  }

  return(f_info)


}
