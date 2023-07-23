test_that("analytical gradient matches numerical gradient when j != jstar", {

  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  z <- rnorm(40) +8
  J <- 10
  p <- 2
  n <- 40
  b0 <- rnorm(10)
  b1 <- 1:10
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = 10, nrow = 40)

  for(i in 1:40){
    for(j in 1:10){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }
  jstar <- 4
  kstar <- 2
  B <- b
  B[1,] <- B[1,] - mean(B[1,])
  B[2,4] <- mean(B[2,-4])
  B[2,] <- B[2,] - mean(B[2,])
  z <- update_z_no_wts(Y,X,B)
  j <- 3
  jstar <- 4
  kstar <- 2

  lstar <- function(Bj){
    B_standin <- B
    B_standin[,j] <- Bj
    B_standin[kstar,jstar] <- mean(B_standin[kstar,-jstar])

    log_mu_j <- X%*%B_standin[,j] + z
    log_mu_jstar <- X%*%B_standin[,jstar] + z
    return( sum(Y[,j]*log_mu_j- exp(log_mu_j)) +
              sum(Y[,jstar]*log_mu_jstar - exp(log_mu_jstar)))
  }

  numerical <- numDeriv::grad(lstar,B[,j])
  dg_dBj <- get_dg_dBj(B,j,mean,kstar)
  fn_analytical <- constr_grad(X,Y,B,z,j,jstar,kstar,dg_dBj)
  expect_equal(as.numeric(fn_analytical),numerical)

  numerical_info <- -1* numDeriv::hessian(lstar,B[,j])

  analytical_info <-constr_info(X,
                                            Y,
                                            B,
                                            z,
                                            j,
                                            jstar,
                                            kstar,
                                            dg_dBj)

  expect_equal(numerical_info,analytical_info)


})

test_that("analytical gradient matches numerical gradient when j = jstar", {

  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  z <- rnorm(40) +8
  J <- 10
  p <- 2
  n <- 40
  b0 <- rnorm(10)
  b1 <- 1:10
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = 10, nrow = 40)

  for(i in 1:40){
    for(j in 1:10){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }
  jstar <- 4
  kstar <- 2
  B <- b
  B[1,] <- B[1,] - mean(B[1,])
  B[2,4] <- mean(B[2,-4])
  B[2,] <- B[2,] - mean(B[2,])
  z <- update_z_no_wts(Y,X,B)
  j <- 4
  jstar <- 4
  kstar <- 2

  lstar <- function(Bj){
    B_standin <- B
    B_standin[,j] <- Bj
    B_standin[kstar,j] <- B[kstar,j]

    log_mu_j <- X%*%B_standin[,j] + z
    return( sum(Y[,j]*log_mu_j- exp(log_mu_j)))
    }

  numerical <- numDeriv::grad(lstar,B[,j])


  dg_dBj <- get_dg_dBj(B,j,mean,kstar)
  fn_analytical <- constr_grad(X,Y,B,z,j,jstar,kstar,dg_dBj)


  expect_equal(as.numeric(fn_analytical),numerical)




})

