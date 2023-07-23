

get_A_n <- function(X_cup,
                    B_cup,
                    Y,
                    J,
                    p,
                    n){

  A <- matrix(0,ncol = p*J,nrow = p*J)
  for(i in 1:n){
    X_cup_i <- X_cup[(i-1)*J + 1:J,]
    z_i <- log(sum(Y[i,])) - log(sum( exp(X_cup%*%B_cup)))
    mu_i <- as.numeric((exp(X_cup_i%*%B_cup + z_i)))
    p_i <- mu_i/sum(mu_i)
    M_i <- outer(rep(1,J),p_i)
    I_minus_M_i <- -1*M_i + diag(J)

    premult <- I_minus_M_i%*%X_cup_i
    A <- A + Matrix::crossprod(premult,diag(mu_i))%*%premult
  }
  return(A)
}
