

micro_nonsimplified_score_stat <- function(Y,
                             X,
                             B_fitted,
                             k,
                             j){
  J <- ncol(Y)
  n <- nrow(Y)
  X_cup <- X_cup_from_X(X,J)
  B_cup <- B_cup_from_B(B_fitted)
  p <- ncol(X)

  UUT <- matrix(0,ncol = p*J, nrow = p*J)
  U <- matrix(0,ncol = p*J, nrow = 1)

  for(i in 1:n){
    X_cup_i <- X_cup[(i-1)*J + 1:J,]
    z_i = log(sum(Y[i,])/sum(exp(X_cup_i%*%B_cup)))
    mu_i = exp(X_cup_i%*%B_cup + z_i)
    p_i = exp(X_cup_i%*%B_cup)
    p_i = p_i/sum(p_i)
    M_i <- -matrix(1,ncol = 1, nrow = J)%*%t(as.matrix(p_i))
    U_i <- t(as.matrix(Y[i,] - mu_i))%*%(diag(rep(1,J)) + M_i)%*%X_cup_i
    U <- U + U_i
    UUT <- UUT + crossprod(as.matrix(U_i))
    # D_i <- diag(as.numeric(mu_i))%*%(diag(rep(1,J)) + M_i)%*%X_cup_i
    # A <- A + Matrix::t(D_i)%*%diag(1/as.numeric(mu_i))%*%D_i
  }

  U <- t(as.matrix(U))

  null_index <- (j - 1)*p + k
  # U[-null_index,] = 0
  return(t(U)%*%MASS::ginv(UUT)%*%U)

}
