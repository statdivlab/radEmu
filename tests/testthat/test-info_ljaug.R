test_that("expected info for augmented ll in Bj equal to numerical hessian when |Bk_constr_j| > delta", {
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
  numeric_hess <- numDeriv::hessian(ljaug,B[,3])
  analytic_info <- info_ljaug(Y = Y,
                              X = X,
                              B = B,
                              z = z,
                              constraint_fn = constraint_fn,
                              constraint_grad_fn = constraint_grad_fn,
                              j = 3, j_constr = 2,k_constr = 2)

  expect_equal( -1* numeric_hess,  analytic_info)
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
  numeric_hess <- numDeriv::hessian(ljaug,B[,19])
  analytic_info <- info_ljaug(Y = Y,
                              X = X,
                              B = B,
                              z = z,
                              constraint_fn = constraint_fn,
                              constraint_grad_fn = constraint_grad_fn,
                              j = 19, j_constr = 2,k_constr = 2)

  expect_equal( -1* numeric_hess,  analytic_info)

})



