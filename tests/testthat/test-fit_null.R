test_that("we get same null fit with different j_ref", {
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
  constraint_fn <- function(x){mean(x)}
  
  ##### Arguments to fix:
  
  # constraint_grad_fn <- function(x){dpseudohuber_center_dx(x,0.1)
  constraint_grad_fn <- function(x){ rep(1/length(x), length(x))}
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
                           verbose = FALSE)#)
  
  B <- full_fit$B
  Y_aug <- full_fit$Y_augmented
  
  
  
  X_cup <- X_cup_from_X(X,J)
  
  fits <- vector(5,mode = "list")
  for(j_ref in 1:5){
    fits[[j_ref]] <- fit_null(B = B,
                              Y = Y_aug,
                              X = X ,
                              X_cup = X_cup,
                              k_constr = k_constr,
                              j_constr = j_constr,
                              j_ref = j_ref,
                              constraint_fn = constraint_fn,
                              constraint_tol = 1e-5,
                              B_tol = 1e-4,
                              constraint_grad_fn = constraint_grad_fn,
                              verbose = FALSE)
  }
  
  Bs <- lapply(fits,function(x){
    temp_B <- x$B
    for(k in 1:2){
      temp_B[k,] <- temp_B[k,] - constraint_fn(temp_B[k,])
    }
    return(B_cup_from_B(temp_B))
  })
  
  expect_equal(Bs[[1]],Bs[[2]],tolerance = 1e-2)
  expect_equal(Bs[[2]],Bs[[3]],tolerance = 1e-2)
  expect_equal(Bs[[3]],Bs[[4]],tolerance = 1e-2)
  expect_equal(Bs[[4]],Bs[[5]],tolerance = 1e-2)
  expect_equal(mean(fits[[1]]$B[2,]), fits[[1]]$B[2,1],tolerance = 1e-4)
})


test_that("Null fit satisfies null constraints", {
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
  
  constraint_fn <- function(x){ pseudohuber_center(x,0.1)}
  # constraint_fn <- function(x){mean(x)}
  
  ##### Arguments to fix:
  
  constraint_grad_fn <- function(x){dpseudohuber_center_dx(x,0.1)}
  # constraint_grad_fn <- function(x){ rep(1/length(x), length(x))}
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
                           tolerance = 1e-5,
                           verbose = FALSE)#)
  
  B <- full_fit$B
  Y_aug <- full_fit$Y_augmented
  
  
  
  X_cup <- X_cup_from_X(X,J)
  
  
  fitted_null <- fit_null(B = B,
                          Y = Y_aug,
                          X = X ,
                          X_cup = X_cup,
                          k_constr = k_constr,
                          j_constr = j_constr,
                          j_ref = 5,
                          constraint_fn = constraint_fn,
                          constraint_tol = 1e-5,
                          B_tol = 1e-1,
                          constraint_grad_fn = constraint_grad_fn,
                          verbose = TRUE)
  
  
  
  expect_true(abs(fitted_null$B[k_constr,j_constr] -
                    constraint_fn(fitted_null$B[k_constr]))<1e-5)
})



test_that("If trackB = TRUE we get B at each iteration back", {
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
  
  constraint_fn <- function(x){ pseudohuber_center(x,0.1)}
  # constraint_fn <- function(x){mean(x)}
  
  ##### Arguments to fix:
  
  constraint_grad_fn <- function(x){dpseudohuber_center_dx(x,0.1)}
  # constraint_grad_fn <- function(x){ rep(1/length(x), length(x))}
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
                           tolerance = 1e-5,
                           verbose = FALSE)#)
  
  B <- full_fit$B
  Y_aug <- full_fit$Y_augmented
  
  
  
  X_cup <- X_cup_from_X(X,J)
  
  
  fitted_null <- fit_null(B = B,
                          Y = Y_aug,
                          X = X ,
                          X_cup = X_cup,
                          k_constr = k_constr,
                          j_constr = j_constr,
                          j_ref = 5,
                          constraint_fn = constraint_fn,
                          constraint_tol = 1e-5,
                          B_tol = 1e-1,
                          constraint_grad_fn = constraint_grad_fn,
                          verbose = FALSE,
                          trackB = TRUE)
  
  expect_true(max(fitted_null$Bs$iter) == 169)
})