test_that("analytical derivative equals numerical derivative", {


  B <- rbind(
    c(0.3838361, -1.512031, -0.3943937, -0.1595407, 0.2784027, 1.3981213, -1.574401, 0.07431514, 1.039499, -0.0320598),
    c(-4.5094179, -3.477205, -2.5046550, -1.5277569, -0.5062945, 0.4977191, 1.571575, 2.49900623, 3.479181, 4.4497043)
  )

  B0 <- rbind(
    c(0.385657, -1.509260, -0.391837, -0.1557335, 0.2858604, 1.4017562, -1.564623, 0.07670448, 1.041326, -0.05164478),
    c(-4.505904, -3.474742, -2.501889, -1.5262399, -0.5084281, 0.4994405, 1.567152, 2.50197755, 3.482716, 4.47465086)
  )

  constraint_fn <- pseudohuber_center
  constraint_grad_fn <- dpseudohuber_center_dx
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
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }
  z <- update_z_no_wts(Y,X,B)
  cfaB0 <- constraint_fn(B0[2,])
  cgaB0 <- constraint_grad_fn(B0[2,])
  llaz <- function(B8){
    temp_B <- B
    temp_B[,8] <- B8
    return(linearized_aug_lag_z(X =X,
                         Y = Y,
                         j = 8,
                         j_constr = 10,
                         k_constr = 2,
                         B = temp_B,
                         z = z,
                         B0 = B0,
                         u = 1e5,
                         rho  = 10,
                         constraint_fn_at_B0 =   cfaB0,
                         constraint_grad_at_B0 =   cgaB0,
                         compute_gradient = FALSE)$value)

  }

  numeric_deriv <- numDeriv::grad(llaz,B[,8])



  analytical_deriv <-
    linearized_aug_lag_z(X =X,
                         Y = Y,
                         j = 8,
                         j_constr = 10,
                         k_constr = 2,
                         B = B,
                         z = z,
                         B0 = B0,
                         u = 1e5,
                         rho  = 10,
                         constraint_fn_at_B0 =   cfaB0,
                         constraint_grad_at_B0 =   cgaB0,
                         compute_gradient = TRUE)

  expect_equal(numeric_deriv,analytical_deriv$gr,tolerance = 0.0001)

  bigger_numeric_deriv <- numDeriv::grad(llaz,B[,8] + c(2,0))

  changed_B <- B
  changed_B[,8] <- changed_B[,8] + c(2,0)
  bigger_analytical_deriv <-
    linearized_aug_lag_z(X =X,
                         Y = Y,
                         j = 8,
                         j_constr = 10,
                         k_constr = 2,
                         B = changed_B,
                         z = z,
                         B0 = B0,
                         u = 1e5,
                         rho  = 10,
                         constraint_fn_at_B0 =   cfaB0,
                         constraint_grad_at_B0 =   cgaB0,
                         compute_gradient = TRUE)


  expect_equal(bigger_numeric_deriv,bigger_analytical_deriv$gr)


})

test_that("analytical hessian equals numerical hessian", {


  B <- rbind(
    c(0.3838361, -1.512031, -0.3943937, -0.1595407, 0.2784027, 1.3981213, -1.574401, 0.07431514, 1.039499, -0.0320598),
    c(-4.5094179, -3.477205, -2.5046550, -1.5277569, -0.5062945, 0.4977191, 1.571575, 2.49900623, 3.479181, 4.4497043)
  )

  B0 <- rbind(
    c(0.385657, -1.509260, -0.391837, -0.1557335, 0.2858604, 1.4017562, -1.564623, 0.07670448, 1.041326, -0.05164478),
    c(-4.505904, -3.474742, -2.501889, -1.5262399, -0.5084281, 0.4994405, 1.567152, 2.50197755, 3.482716, 4.47465086)
  )

  constraint_fn <- pseudohuber_center
  constraint_grad_fn <- dpseudohuber_center_dx
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
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }
  z <- update_z_no_wts(Y,X,B)
  cfaB0 <- constraint_fn(B0[2,])
  cgaB0 <- constraint_grad_fn(B0[2,])
  llaz <- function(B8){
    temp_B <- B
    temp_B[,8] <- B8
    return(linearized_aug_lag_z(X =X,
                                Y = Y,
                                j = 8,
                                j_constr = 10,
                                k_constr = 2,
                                B = temp_B,
                                z = z,
                                B0 = B0,
                                u = 0,
                                rho  = 0,
                                constraint_fn_at_B0 =   cfaB0,
                                constraint_grad_at_B0 =   cgaB0,
                                compute_gradient = FALSE)$value)

  }

  numeric_hess <- numDeriv::hessian(llaz,B[,8])



  analytical_hess <-
    linearized_aug_lag_z(X =X,
                         Y = Y,
                         j = 8,
                         j_constr = 10,
                         k_constr = 2,
                         B = B,
                         z = z,
                         B0 = B0,
                         u = 0,
                         rho  = 0,
                         constraint_fn_at_B0 =   cfaB0,
                         constraint_grad_at_B0 =   cgaB0,
                         compute_hessian = TRUE)

  expect_equal(numeric_hess,analytical_hess$hess)

  bigger_numeric_hess <- numDeriv::hessian(llaz,B[,8] + c(2,0))

  changed_B <- B
  changed_B[,8] <- changed_B[,8] + c(2,0)
  bigger_analytical_hess <-
    linearized_aug_lag_z(X =X,
                         Y = Y,
                         j = 8,
                         j_constr = 10,
                         k_constr = 2,
                         B = changed_B,
                         z = z,
                         B0 = B0,
                         u = 1e5,
                         rho  = 10,
                         constraint_fn_at_B0 =   cfaB0,
                         constraint_grad_at_B0 =   cgaB0,
                         compute_hessian = TRUE)


  expect_equal(bigger_numeric_hess,bigger_analytical_hess$hess)


})
