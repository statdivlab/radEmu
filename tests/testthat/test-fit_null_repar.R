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
