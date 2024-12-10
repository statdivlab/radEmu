test_that("we get same null fit using fit_null and fit_null_repar", {
  set.seed(59542234)
  n <- 10
  p <- 2
  X <- cbind(1,rep(c(0,1),each = n/2))
  J <- 5
  z <- rnorm(n) +8
  b0 <- rnorm(J)
  b1 <- seq(1,5,length.out = J)
  b1 <- b1 - mean(b1)
  b0 <- b0 - mean(b0)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)
  
  k_constr <- 2
  j_constr <- 1
  p <- 2
  
  # constraint_fn <- function(x){ pseudohuber_center(x,0.1)}
  constraint_fn <- function(x){pseudohuber_center(x,0.1)}
  
  ##### Arguments to fix:
  
  constraint_grad_fn <- function(x){dpseudohuber_center_dx(x,0.1)}

  rho_init = 1
  tau = 1.2
  kappa = 0.8
  obj_tol = 100
  B_tol <- 1e-3
  constraint_tol = 1e-5
  init_tol = 1e6
  c1 = 1e-4
  maxit = 1000
  inner_maxit = 25
  
  Y[] <- 0
  for(i in 1:n){
    while(sum(Y[i,])==0){
      for(j in 1:J){
        temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
        Y[i,j] <- rpois(1, lambda = temp_mean)
        # Y[i,j] <- rnbinom(1,mu = temp_mean, size = 3)*rbinom(1,1,0.6)
      }
    }
  }
  
  full_fit <- #suppressMessages(
    emuFit_micro_penalized(X = X,
                           Y = Y,
                           B = NULL,
                           constraint_fn = mean,
                           tolerance = 1e-7,
                           verbose = FALSE,
                           j_ref = 5)
  
  B <- full_fit$B
  Y_aug <- full_fit$Y_augmented
  
  X_cup <- X_cup_from_X(X,J)
 
  null_fit <-  fit_null(B = B,
                        Y = Y_aug,
                        X = X ,
                        X_cup = X_cup,
                        k_constr = k_constr,
                        j_constr = j_constr,
                        j_ref = 5,
                        constraint_fn = constraint_fn,
                        constraint_tol = 1e-4,
                        B_tol = 1e-4,
                        constraint_grad_fn = constraint_grad_fn,
                        verbose = FALSE,
                        trackB = FALSE)
  
  null_fit_repar <-
    fit_null_repar(B = B,
                   Y = Y_aug,
                   X = X,
                   k_constr = k_constr,
                   j_constr = j_constr,
                   j_ref =5,
                   constraint_fn = constraint_fn,
                   constraint_grad_fn = constraint_grad_fn,
                   maxit = 1000,
                   verbose = FALSE,
                   trackB = TRUE,
                   tolerance = 1e-2,
                   method = "fisher",
                   starting_stepsize = 0.5)
  
  expect_equal(null_fit_repar$B,null_fit$B,tolerance = 0.1)
})




test_that("Vast majority of deviation from zero in null fits computed via 
fit_null and fit_null_repar (with continuous covariate) are explained by 
derivative of constraint function (as they should be at a constrained 
approximate optimum).", {
  set.seed(59542234)
  n <- 10
  p <- 2
  X <- cbind(1,rnorm(10))
  J <- 25
  z <- rnorm(n) +8
  b0 <- rnorm(J)
  b1 <- seq(1,5,length.out = J)
  b1 <- b1 - mean(b1)
  # b1[1] <- pseudohuber_center(b1[-1],0.1)
  b0 <- b0 - mean(b0)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)

  k_constr <- 2
  j_constr <- 1
  p <- 2

  constraint_fn <- function(x){ pseudohuber_center(x,0.1)}
  # constraint_fn <- function(x){pseudohuber_center(x,0.1)}
  # constraint_fn <- function(x){mean(x)}

  ##### Arguments to fix:

  constraint_grad_fn <- function(x){dpseudohuber_center_dx(x,0.1)}
  # constraint_grad_fn <- function(x){rep(1/length(x),length(x))}

  rho_init = 1
  tau = 1.2
  kappa = 0.8
  obj_tol = 100
  B_tol <- 1e-3
  constraint_tol = 1e-5
  init_tol = 1e6
  c1 = 1e-4
  maxit = 1000
  inner_maxit = 25

  Y[] <- 0
  for(i in 1:n){
    while(sum(Y[i,])==0){
      for(j in 1:J){
        temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
        Y[i,j] <- rpois(1, lambda = temp_mean)
        # Y[i,j] <- rnbinom(1,mu = temp_mean, size = 3)*rbinom(1,1,0.6)
      }
    }
  }

  full_fit <- #suppressMessages(
    emuFit_micro_penalized(X = X,
                           Y = Y,
                           B = NULL,
                           constraint_fn = mean,
                           tolerance = 1e-6,
                           verbose = FALSE,
                           j_ref = 25)

  B <- full_fit$B
  Y_aug <- full_fit$Y_augmented

  X_cup <- X_cup_from_X(X,J)

  null_fit <-  fit_null(B = B,
                        Y = Y_aug,
                        X = X ,
                        X_cup = X_cup,
                        k_constr = k_constr,
                        j_constr = j_constr,
                        j_ref = 25,
                        constraint_fn = constraint_fn,
                        constraint_tol = 1e-3,
                        B_tol = 1e-2,
                        rho_init = 100,
                        kappa = kappa,
                        tau = 1.2,
                        inner_tol = 0.01,
                        inner_maxit = 50,
                        maxit = 1e4,
                        constraint_grad_fn = constraint_grad_fn,
                        verbose = FALSE,
                        trackB = TRUE)

  null_fit_repar <-
    fit_null_repar(B = B,
                   Y = Y_aug,
                   X = X,
                   k_constr = k_constr,
                   j_constr = j_constr,
                   j_ref =25,
                   constraint_fn = constraint_fn,
                   constraint_grad_fn = constraint_grad_fn,
                   maxit = 10000,
                   verbose = FALSE,
                   trackB = TRUE,
                   tolerance = 1e-4,
                   method = "fisher",
                   starting_stepsize = 1,
                   do_shift = TRUE)


  
  null_fit_grad <- numDeriv::grad(function(x){
    profile_ll(X = X, Y = Y_aug, B = cbind(rbind(x[1:(J - 1)],
                                                 x[(J -1) + 1:(J - 1)]),
                                           0))
  },
  c(null_fit$B[1,1:(J - 1)],null_fit$B[2,1:(J - 1)])
  )
  
  null_fit_repar_grad <- numDeriv::grad(function(x){
    profile_ll(X = X, Y = Y_aug, B = cbind(rbind(x[1:(J - 1)],
                                                 x[(J -1) + 1:(J - 1)]),
                                           0))
  },
  c(null_fit_repar$B[1,1:(J - 1)],null_fit_repar$B[2,1:(J - 1)])
  )
  
  cg_term <- constraint_grad_fn(null_fit$B[2,])
  cg_term[1] <- cg_term[1] - 1
  cg_term <- cg_term[1:24]
  nf_b2 <- null_fit_grad[25:48]
  cg_lm <- lm( nf_b2~cg_term - 1)
 

  
  cg_repar_term <- constraint_grad_fn(null_fit_repar$B[2,])
  cg_repar_term[1] <- cg_repar_term[1] - 1
  cg_repar_term <- cg_repar_term[1:24]
  
  nfr_b2 <- null_fit_repar_grad[25:48]
  cg_repar_lm <- lm( nfr_b2~cg_repar_term - 1)

  
  expect_equal(summary(cg_lm)$r.squared,1,tolerance  = 0.0001)
  expect_equal(summary(cg_repar_lm)$r.squared,1,tolerance  = 0.0001)
})
