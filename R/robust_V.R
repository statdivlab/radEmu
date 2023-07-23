
robust_V <- function(X,Y,B,M = 50){
  J <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  X_cup <- X_cup_from_X(X,J)
  B_cup <- B_cup_from_B(B)
  Uis <- matrix(0,nrow = p*J,ncol = n)

  for(i in 1:n){
    X_cup_i <- X_cup[(i-1)*J + 1:J,]
    z_i = log(sum(Y[i,])/sum(exp(X_cup_i%*%B_cup)))
    mu_i = exp(X_cup_i%*%B_cup + z_i)
    p_i = exp(X_cup_i%*%B_cup)
    p_i = p_i/sum(p_i)
    M_i <- -matrix(1,ncol = 1, nrow = J)%*%t(as.matrix(p_i))
    Uis[,i] <- as.numeric(t(as.matrix(Y[i,] - mu_i))%*%(diag(rep(1,J)) + M_i)%*%X_cup_i)
  }

  V_est <- matrix(0,ncol = p*J,nrow = p*J)

  for(m in 1:M){
    group1 <- sample(1:n,ceiling(n/2),replace = FALSE)

    UUT1 <- UUT2 <- matrix(0,nrow = p*J, ncol = p*J)

    for(i in 1:n){
      if(i %in% group1){
        UUT1 <- UUT1 + tcrossprod(Uis[,i,drop = FALSE])
      } else{
        UUT2 <- UUT2 + tcrossprod(Uis[,i,drop = FALSE])
      }
    }
    UUT1 <- UUT1/length(group1)
    UUT2 <- UUT2/(n - length(group1))
    UUT1_eigen <- eigen(UUT1)
    V_est <- V_est + UUT1_eigen$vectors%*%(t(UUT1_eigen$vectors)%*%UUT2%*%UUT1_eigen$vectors)%*%t(UUT1_eigen$vectors)
  }
  V_est <- (1/M)*V_est

  return(V_est)
}
