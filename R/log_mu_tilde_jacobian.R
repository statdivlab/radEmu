

log_mu_tilde_jacobian <- function(X_tilde,mu_tilde,n,J){

jacobian <- 0*X_tilde

for(i in 1:n){
  print(i)
  A_mu <- matrix(rep(mu_tilde[(i - 1)*J + 1:J],J),ncol = J,byrow = TRUE)
  jacobian[(i - 1)*J + 1:J,] <- X_tilde[(i - 1)*J + 1:J,] -
    A_mu%*%X_tilde[(i - 1)*J + 1:J,]
}
}
