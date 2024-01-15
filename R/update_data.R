

update_data <- function(Y,
                        X_tilde,
                        B,
                        p,
                        n,
                        J,
                        verbose = FALSE){
  X_tilde_J <- X_tilde
  indexes_to_zero <- numeric(0)
  for(i in 1:n){
    indexes_to_zero <- c(indexes_to_zero,(i - 1)*J + J)
  }
  X_tilde_J[indexes_to_zero,] <- 0
  # for(i in 1:n){
  #   X_tilde_J[(i - 1)*J + J,] <- 0#*X_tilde_J[(i - 1)*J + J,]
  # }
  X_tilde_J <- X_tilde_J[,1:(p*(J - 1))]
  B_J <- B
  for(k in 1:p){
    B_J[k,] <- B_J[k,] - B_J[k,J]
  }
  beta_tilde_J <- B_cup_from_B(B_J[,1:(J - 1),drop= FALSE])

  tau <- rowSums(Y)
  if(sum(tau ==0)>0){
  tau[tau == 0] <- 0.1
  warning("At least one observation has total count equal to zero!")
  }
    # sapply(1:n,
    #        function(i)
    #          log(sum(Y[i,])) -
    #          log(sum(exp(X_tilde_J[(i - 1)*J + 1:J,]%*%beta_tilde_J))))



  augmentations <-
    asymmetric_hat(beta_tilde_J = beta_tilde_J,
                   tau = tau,
                   X_tilde_J = X_tilde_J,
                   J = J,
                   n = n,
                   p = p,
                   verbose = verbose)

  return(Y + augmentations)
}

