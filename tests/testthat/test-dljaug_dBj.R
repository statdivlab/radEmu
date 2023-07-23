test_that("derivative of augmented ll in Bj correct when |Bk_constr_j| > delta", {
  set.seed(4323)
  J = 40
  p = 2
  n <- 100
  X <- cbind(1,rep(c(0,1),each = n/2))
  z <- rnorm(n)
  # b0 <- c(-1,0,1)
  # b1 <- c(-1,0,1)
  b0 <- rnorm(J)
  b1 <- seq(-3,3,length.out = J)
  # b1[2] <- mean(b1[-2])
  b1[2] <- huber_center(b1[-2])
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)


  for(i in 1:n){
    for(j in 1:J){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
      # Y[i,j] <- rnbinom(1,size = 5,mu = temp_mean)
    }
  }

  B <- rbind(b0,b1)

  z <- update_z_no_wts(Y,X,B)

  constraint_fn <- huber_center

  constraint_grad_fn <- function(x){
    hc <- huber_center(x)
    unscaled_deriv <- as.numeric(abs(x - hc)<=1)
    return(unscaled_deriv/sum(unscaled_deriv))
  }

  numeric_constraint_deriv <- numDeriv::grad(constraint_fn, B[2,])
  analytic_constraint_deriv <- constraint_grad_fn(B[2,])
  expect_true(sum((numeric_constraint_deriv - analytic_constraint_deriv)^2)<1e-12)


  ljaug <- function(Bj){
    lmj <- X%*%matrix(Bj,ncol = 1) + z
    llj <- sum(Y[,j,drop = FALSE]*lmj - exp(lmj))

    Bj_constr <- B[,j_constr]
    Bk_constr <- B[k_constr,]
    Bk_constr[j] <- Bj[k_constr]
    Bj_constr[k_constr] <- constraint_fn(Bk_constr[-j_constr])

    lmj_constr <- X%*%matrix(Bj_constr,ncol = 1) + z
    llj_constr <- sum(Y[,j_constr,drop = FALSE]*lmj_constr - exp(lmj_constr))

    return(llj + llj_constr)

  }

  j <- 3
  j_constr <- 2
  k_constr <- 2
  numeric_grad <- numDeriv::grad(ljaug,B[,3])
  analytic_grad <- dljaug_dBj(Y = Y,
                              X = X,
                              B = B,
                              z = z,
                              constraint_fn = constraint_fn,
                              constraint_grad_fn = constraint_grad_fn,
                              j = 3, j_constr = 2,k_constr = 2)

 expect_equal(numeric_grad,as.numeric(analytic_grad))
})


test_that("derivative of augmented ll in Bj correct when |Bk_constr_j| <= delta", {
  set.seed(4323)
  J = 40
  p = 2
  n <- 100
  X <- cbind(1,rep(c(0,1),each = n/2))
  z <- rnorm(n)
  # b0 <- c(-1,0,1)
  # b1 <- c(-1,0,1)
  b0 <- rnorm(J)
  b1 <- seq(-3,3,length.out = J)
  # b1[2] <- mean(b1[-2])
  b1[2] <- huber_center(b1[-2])
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)


  for(i in 1:n){
    for(j in 1:J){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
      # Y[i,j] <- rnbinom(1,size = 5,mu = temp_mean)
    }
  }

  B <- rbind(b0,b1)

  z <- update_z_no_wts(Y,X,B)

  constraint_fn <- huber_center

  constraint_grad_fn <- function(x){
    hc <- huber_center(x)
    unscaled_deriv <- as.numeric(abs(x - hc)<=1)
    return(unscaled_deriv/sum(unscaled_deriv))
  }

  numeric_constraint_deriv <- numDeriv::grad(constraint_fn, B[2,])
  analytic_constraint_deriv <- constraint_grad_fn(B[2,])
  expect_true(sum((numeric_constraint_deriv - analytic_constraint_deriv)^2)<1e-12)


  ljaug <- function(Bj){
    lmj <- X%*%matrix(Bj,ncol = 1) + z
    llj <- sum(Y[,j,drop = FALSE]*lmj - exp(lmj))

    Bj_constr <- B[,j_constr]
    Bk_constr <- B[k_constr,]
    Bk_constr[j] <- Bj[k_constr]
    Bj_constr[k_constr] <- constraint_fn(Bk_constr[-j_constr])

    lmj_constr <- X%*%matrix(Bj_constr,ncol = 1) + z
    llj_constr <- sum(Y[,j_constr,drop = FALSE]*lmj_constr - exp(lmj_constr))

    return(llj + llj_constr)

  }

  j <- 19
  j_constr <- 2
  k_constr <- 2
  numeric_grad <- numDeriv::grad(ljaug,B[,19])
  analytic_grad <- dljaug_dBj(Y = Y,
                              X = X,
                              B = B,
                              z = z,
                              constraint_fn = constraint_fn,
                              constraint_grad_fn = constraint_grad_fn,
                              j = 19, j_constr = 2,k_constr = 2)

  expect_equal(numeric_grad,as.numeric(analytic_grad))
})


test_that("derivative is right in sitauation we've been stalling in",{

  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  z <- rnorm(40) +5
  J <- 10
  p <- 2
  n <- 40
  b0 <- rnorm(J)
  b1 <- seq(1,10,length.out = J)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = 40)

  for(i in 1:40){
    for(j in 1:J){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rnbinom(1, mu = temp_mean,size = 1.2)#rpois(1, lambda = temp_mean)
    }
  }


  B <- rbind(
    c(-24.08636,-1.259031,-0.1416076,0.09449586,0.5360898,1.65198559,-1.314394,0.3269339,1.291555,1.985846e-01),
  c(19.44656,-3.998268,-3.0254151,-2.04976587,-1.0319541,-0.02408555,1.043626,1.9784515,2.959190,-2.484989e-09)
  )

z <- c(2.177852,4.044364,3.643078,4.863511,5.082759,5.750786,5.315855,5.241445,4.723171,6.621183,4.927388,3.575811,6.803497,3.772744,4.649128,4.500806,
4.452246,5.656701,4.516250,4.234010,10.992182,10.121600,10.822365,10.148005,12.743986,10.734118,12.221640,11.751059,11.890335,11.882248,10.646250,11.987260,
9.817554,12.873551,9.467366,10.592537,12.328858,12.726294,11.779928,9.121726)
constraint_fn <- pseudohuber_center
constraint_grad_fn <- dpseudohuber_center_dx
j <- 1
j_constr <- 10
k_constr <- 2
analytic_grad <- dljaug_dBj(Y = Y,
                            X = X,
                            B = B,
                            z = z,
                            constraint_fn = constraint_fn,
                            constraint_grad_fn = constraint_grad_fn,
                            j = 1, j_constr = 10,k_constr = 2)



ljaug <- function(Bj){
  lmj <- X%*%matrix(Bj,ncol = 1) + z
  llj <- sum(Y[,j,drop = FALSE]*lmj - exp(lmj))

  Bj_constr <- B[,j_constr]
  Bk_constr <- B[k_constr,]
  Bk_constr[j] <- Bj[k_constr]
  Bj_constr[k_constr] <- constraint_fn(Bk_constr[-j_constr])

  lmj_constr <- X%*%matrix(Bj_constr,ncol = 1) + z
  llj_constr <- sum(Y[,j_constr,drop = FALSE]*lmj_constr - exp(lmj_constr))

  return(llj + llj_constr)

}


ll_fn <- function(Bj){
  temp_B <- B
  temp_B[,j] <- Bj
  temp_B[k_constr,j_constr] <- constraint_fn(temp_B[k_constr,-j_constr])
  temp_means <- exp(X%*%temp_B + matrix(z,ncol = 1)%*%matrix(1,nrow = 1, ncol = J))
  return(sum(Y*log(temp_means) - temp_means))
}

numerical_grad <- numDeriv::grad(ljaug,B[,1])
numerical_grad_2 <- numDeriv::grad(ll_fn,B[,1])

expect_true(max(abs(analytic_grad - numerical_grad))<1e-2)
expect_true(max(abs(analytic_grad - numerical_grad_2))<1e-2)

})





