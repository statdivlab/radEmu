test_that("We take same step as we'd take using numerical derivatives when gap, rho, u are zero", {
  set.seed(59542234)
  n <- 10
  p <- 2
  X <- cbind(1,rep(c(0,1),each = n/2))
  J <- 5
  z <- rnorm(n) +8
  b0 <- rnorm(J)
  b1 <- seq(1,10,length.out = J)
  b1 <- b1 - mean(b1)
  b0 <- b0 - mean(b0)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)

  k_constr <- 2
  j_constr <- 1
  p <- 2

  constraint_fn <- rep(list(function(x){ pseudohuber_median(x,0.1)}), 2)
  # constraint_fn <- function(x){mean(x)}

  ##### Arguments to fix:

  constraint_grad_fn <- rep(list(function(x){dpseudohuber_median_dx(x,0.1)}), 2)
  # constraint_grad_fn <- function(x){ rep(1/length(x), length(x))}
  rho_init = 1
  tau = 1.2
  kappa = 0.8
  obj_tol = 100
  score_tol <- 1e-3
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

  j_ref <- 5

  full_fit <- #suppressMessages(
    emuFit_micro_penalized(X = X,
                           Y = Y,
                           B = NULL,
                           constraint_fn = rep(list(mean), 2), 
                           tolerance = 1e-5,
                           verbose = FALSE)#)

  B <- full_fit$B
  B[1,] <- B[1,] - B[1,j_ref]
  B[2,] <- B[2,] - B[2,j_ref]

  Y_aug <- full_fit$Y_augmented

  B[k_constr,j_constr] <- constraint_fn[[k_constr]](B[k_constr,-j_constr])

  u <- 0
  rho <- 0

  z <- update_z(Y,X,B)


  update <- macro_fisher_null(X = X,
                    Y = Y_aug,
                    B = B,
                    z = z,
                    J = J,
                    p = p,
                    k_constr = k_constr,
                    j_constr = j_constr,
                    j_ref = j_ref,
                    rho = rho,
                    u = u,
                    constraint_fn = constraint_fn,
                    constraint_grad_fn = constraint_grad_fn,
                    stepsize = 0.5,
                    c1 = 1e-4,
                    regularization = 0,
                    debug= TRUE)

  lag_fn <- function(b){
    tempB <- B_from_B_cup(b,J,p)
    tempB[,j_ref] <- 0
    gap <- constraint_fn[[k_constr]](tempB[k_constr,]) - tempB[k_constr,j_constr]
    log_means <- do.call(cbind,
                         lapply(1:J,
                                function(j) z + X%*%tempB[,j]))
    return(0.5*rho*gap^2 + u*gap -sum(Y_aug*log_means - exp(log_means)))
  }

  nd <- numDeriv::grad(lag_fn,x= B_cup_from_B(B))
  nh <- numDeriv::hessian(lag_fn,x= B_cup_from_B(B))
  nd <- nd[-((j_ref- 1)*p + 1:p)]
  nh <- nh[-((j_ref- 1)*p + 1:p),-((j_ref- 1)*p + 1:p)]

  n_update <- B_from_B_cup(qr.solve(nh,nd), J = J, p= p)
  max_ratio <- max(update$update/n_update,na.rm = TRUE)

  min_ratio <- min(update$update/n_update,na.rm = TRUE)

  expect_equal(max_ratio,min_ratio,tolerance = 1e-2)

})


test_that("We take same step as we'd take using numerical derivatives when gap is zero and u is zero", {
  set.seed(5234)
  n <- 10
  p <- 2
  X <- cbind(1,rep(c(0,1),each = n/2))
  J <- 5
  z <- rnorm(n) +8
  b0 <- rnorm(J)
  b1 <- seq(1,10,length.out = J)
  b1 <- b1 - mean(b1)
  b0 <- b0 - mean(b0)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)

  k_constr <- 2
  j_constr <- 1
  p <- 2

  constraint_fn <- rep(list(function(x){ pseudohuber_median(x,0.1)}), 2)
  # constraint_fn <- function(x){mean(x)}

  ##### Arguments to fix:

  constraint_grad_fn <- rep(list(function(x){dpseudohuber_median_dx(x,0.1)}), 2)
  # constraint_grad_fn <- function(x){ rep(1/length(x), length(x))}
  rho_init = 1
  tau = 1.2
  kappa = 0.8
  obj_tol = 100
  score_tol <- 1e-3
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

  j_ref <- 5

  full_fit <- #suppressMessages(
    emuFit_micro_penalized(X = X,
                           Y = Y,
                           B = NULL,
                           constraint_fn = rep(list(mean), 2), 
                           tolerance = 1e-5,
                           verbose = FALSE)#)

  B <- full_fit$B
  B[1,] <- B[1,] - B[1,j_ref]
  B[2,] <- B[2,] - B[2,j_ref]

  Y_aug <- full_fit$Y_augmented

  B[k_constr,j_constr] <- constraint_fn[[k_constr]](B[k_constr,-j_constr])

  u <- 0
  rho <- 100

  z <- update_z(Y,X,B)


  update <- macro_fisher_null(X = X,
                              Y = Y_aug,
                              B = B,
                              z = z,
                              J = J,
                              p = p,
                              k_constr = k_constr,
                              j_constr = j_constr,
                              j_ref = j_ref,
                              rho = rho,
                              u = u,
                              constraint_fn = constraint_fn,
                              constraint_grad_fn = constraint_grad_fn,
                              stepsize = 0.5,
                              c1 = 1e-4,
                              regularization = 0,
                              debug= TRUE)

  lag_fn <- function(b){
    tempB <- B_from_B_cup(b,J,p)
    tempB[,j_ref] <- 0
    gap <- constraint_fn[[k_constr]](tempB[k_constr,]) - tempB[k_constr,j_constr]
    log_means <- do.call(cbind,
                         lapply(1:J,
                                function(j) z + X%*%tempB[,j]))
    return(0.5*rho*gap^2 + u*gap -sum(Y_aug*log_means - exp(log_means)))
  }

  nd <- numDeriv::grad(lag_fn,x= B_cup_from_B(B))
  nh <- numDeriv::hessian(lag_fn,x= B_cup_from_B(B))
  nd <- nd[-((j_ref- 1)*p + 1:p)]
  nh <- nh[-((j_ref- 1)*p + 1:p),-((j_ref- 1)*p + 1:p)]

  expect_equal(nd,as.numeric(as.matrix(update$lag_deriv)))


  expect_equal(nh,solve(as.matrix(update$info_inverse)),
               tolerance = 1e-3)

  n_update <- -0.5*B_from_B_cup(qr.solve(nh,nd), J = J, p= p)
  n_update[is.na(n_update)] <- 0

  expect_equal(update$update,n_update,tolerance = 1e-1)

})



test_that("We take similar step as we'd take using numerical derivatives in a more general setting", {
  set.seed(593234)
  n <- 10
  p <- 2
  X <- cbind(1,rep(c(0,1),each = n/2))
  J <- 5
  z <- rnorm(n) +8
  b0 <- rnorm(J)
  b1 <- seq(1,10,length.out = J)
  b1 <- b1 - mean(b1)
  b0 <- b0 - mean(b0)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)

  k_constr <- 2
  j_constr <- 1
  p <- 2

  # constraint_fn <- function(x){ pseudohuber_median(x,0.1)}
  constraint_fn <- rep(list(function(x){mean(x)}), 2)

  ##### Arguments to fix:

  # constraint_grad_fn <- function(x){dpseudohuber_median_dx(x,0.1)}
  constraint_grad_fn <- rep(list(function(x){ rep(1/length(x), length(x))}), 2)
  rho_init = 1
  tau = 1.2
  kappa = 0.8
  obj_tol = 100
  score_tol <- 1e-3
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

  j_ref <- 5

  full_fit <- #suppressMessages(
    emuFit_micro_penalized(X = X,
                           Y = Y,
                           B = NULL,
                           constraint_fn = rep(list(mean), 2), 
                           tolerance = 1e-5,
                           verbose = FALSE)#)

  B <- full_fit$B
  B[1,] <- B[1,] - B[1,j_ref]
  B[2,] <- B[2,] - B[2,j_ref]

  Y_aug <- full_fit$Y_augmented

  u <- -595999
  rho <- 100000

  z <- update_z(Y,X,B)


  update <- macro_fisher_null(X = X,
                              Y = Y_aug,
                              B = B,
                              z = z,
                              J = J,
                              p = p,
                              k_constr = k_constr,
                              j_constr = j_constr,
                              j_ref = j_ref,
                              rho = rho,
                              u = u,
                              constraint_fn = constraint_fn,
                              constraint_grad_fn = constraint_grad_fn,
                              stepsize = 0.5,
                              c1 = 1e-4,
                              regularization = 0,
                              max_step = 10,
                              debug= TRUE)

  lag_fn <- function(b){
    tempB <- B_from_B_cup(b,J,p)
    tempB[,j_ref] <- 0
    gap <- constraint_fn[[k_constr]](tempB[k_constr,]) - tempB[k_constr,j_constr]
    log_means <- do.call(cbind,
                         lapply(1:J,
                                function(j) z + X%*%tempB[,j]))
    return(0.5*rho*gap^2 + u*gap -sum(Y_aug*log_means - exp(log_means)))
  }

  nd <- numDeriv::grad(lag_fn,x= B_cup_from_B(B))
  nh <- numDeriv::hessian(lag_fn,x= B_cup_from_B(B))
  nd <- nd[-((j_ref- 1)*p + 1:p)]
  nh <- nh[-((j_ref- 1)*p + 1:p),-((j_ref- 1)*p + 1:p)]

  # plot(asinh(nd),asinh(as.numeric(as.matrix(update$lag_deriv))))
  # abline(a = 0, b= 1,lty = 2)

  expect_equal(nd,as.numeric(as.matrix(update$lag_deriv)))


  expect_equal(nh,solve(as.matrix(update$info_inverse)),
               tolerance = 1e-5)

  # plot(asinh(diag(nh)),
  #      asinh(diag(solve(as.matrix(update$info_inverse)))))

  n_update <- -update$stepsize*B_from_B_cup(qr.solve(nh,nd), J = J, p= p)
  n_update[is.na(n_update)] <- 0

  expect_equal(update$update,n_update,tolerance = 1e-1)

})

test_that("We take same step as we'd take using numerical derivatives when gap, rho, u are zero", {
  set.seed(59542234)
  n <- 10
  p <- 2
  X <- cbind(1,rep(c(0,1),each = n/2))
  J <- 10000
  z <- rnorm(n) +8
  b0 <- rnorm(J)
  b1 <- seq(1,10,length.out = J)
  b1 <- b1 - mean(b1)
  b0 <- b0 - mean(b0)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)
  
  k_constr <- 2
  j_constr <- 1
  p <- 2
  
  constraint_fn <- rep(list(function(x){ pseudohuber_median(x,0.1)}), 2)
  # constraint_fn <- function(x){mean(x)}
  
  ##### Arguments to fix:
  
  constraint_grad_fn <- rep(list(function(x){dpseudohuber_median_dx(x,0.1)}), 2)
  # constraint_grad_fn <- function(x){ rep(1/length(x), length(x))}
  rho_init = 1
  tau = 1.2
  kappa = 0.8
  obj_tol = 100
  score_tol <- 1e-3
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
  
  j_ref <- 5
  
  full_fit <- #suppressMessages(
    emuFit_micro_penalized(X = X,
                           Y = Y,
                           B = NULL,
                           constraint_fn = rep(list(mean), 2), 
                           tolerance = 10,
                           verbose = FALSE)#)
  
  B <- full_fit$B
  B[1,] <- B[1,] - B[1,j_ref]
  B[2,] <- B[2,] - B[2,j_ref]
  
  Y_aug <- full_fit$Y_augmented
  
  B[k_constr,j_constr] <- constraint_fn[[k_constr]](B[k_constr,-j_constr])
  
  u <- 0
  rho <- 0
  
  z <- update_z(Y,X,B)
  
  macro_time <- try(system.time(macro_fisher_null(X = X,
                              Y = Y_aug,
                              B = B,
                              z = z,
                              J = J,
                              p = p,
                              k_constr = k_constr,
                              j_constr = j_constr,
                              j_ref = j_ref,
                              rho = rho,
                              u = u,
                              constraint_fn = constraint_fn,
                              constraint_grad_fn = constraint_grad_fn,
                              stepsize = 0.5,
                              c1 = 1e-4,
                              regularization = 0,
                              debug= FALSE)))

  if(inherits(macro_time, "try-error")){
    expect_true(FALSE)
  }
  expect_true(macro_time[3]<15)
  #failing this test indicates that macro_fisher_null may be directly 
  #computing the full inverse of the (approximate) hessian matrix of the log likelihood
  #this includes an outer product that does not need to be computed
})

