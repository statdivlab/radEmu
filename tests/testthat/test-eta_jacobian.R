test_that("multiplication works", {

  Y <- cbind(c(4,6),c(2,3))
  X <- cbind(1,c(0,1))
  X_tilde <- X_cup_from_X(X,2)
  n <- 2
  J <- 2
  p <- 2

  X_tilde_J <- X_tilde
  for(i in 1:n){
    X_tilde_J[(i - 1)*J + J,] <- 0*X_tilde_J[(i - 1)*J + J,]
  }
  X_tilde_J <- X_tilde_J[,1:(p*(J - 1))]

  B <- emuFit_micro(X,Y,constraint_fn = function(x) x[2],tolerance = 1e-5)
  beta_tilde_J <- B[,1,drop = FALSE]
  tau <- rowSums(Y)

  get_eta <- function(beta_tilde_J, tau){
    mu <- exp(X_tilde_J%*%matrix(beta_tilde_J,ncol = 1))
    for(i in 1:n){
      mu[(i - 1)*J + 1:J,] <- mu[(i - 1)*J + 1:J,]/sum(mu[(i - 1)*J + 1:J,])
      mu[(i - 1)*J + 1:J,] <- mu[(i - 1)*J + 1:J,]*tau[i]
    }
    return(as.numeric(as.matrix(log(mu))))
  }




  numerical_deriv_1 <- numDeriv::grad(function(x) get_eta(x[1:2],x[3:4])[1],x = c(as.numeric(beta_tilde_J),tau))
  numerical_deriv_2 <- numDeriv::grad(function(x) get_eta(x[1:2],x[3:4])[2],x = c(as.numeric(beta_tilde_J),tau))

  analytical_deriv <- eta_jacobian(beta_tilde_J = beta_tilde_J,
                           tau = tau,
                           X_tilde_J = X_tilde_J,
                           J = 2,
                           i =1,
                           n = 2)

  numerical_deriv <- rbind(numerical_deriv_1,numerical_deriv_2)
  rownames(numerical_deriv) <- NULL

  expect_equal(numerical_deriv,
               as.matrix(analytical_deriv))


})
