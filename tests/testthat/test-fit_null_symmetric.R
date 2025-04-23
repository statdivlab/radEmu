test_that("we get same null fit with different j_ref", {
  set.seed(59542234)
  n <- 10
  J <- 5
  X <- cbind(1,rep(c(0,1),each = n/2))
  b0 <- rnorm(J)
  b1 <- seq(1,10,length.out = J)
  b1 <- b1 - mean(b1)
  b0 <- b0 - mean(b0)
  Y <- radEmu:::simulate_data(n = n,
                              J = J,
                              X = X,
                              b0 = b0,
                              b1 = b1,
                              distn = "Poisson",
                              mean_z = 8)
  
  k_constr <- 2
  j_constr <- 1
  p <- 2
  
  # constraint_fn <- rep(list(function(x){mean(x)}), 2)
  constraint_fn <- rep(list(function(x)pseudohuber_center(x,0.1)), 2)
  ##### Arguments to fix:
  
  # constraint_grad_fn <- function(x){dpseudohuber_center_dx(x,0.1)
  constraint_grad_fn <- rep(list(function(x){dpseudohuber_center_dx(x,0.1)}), 2)
  
  maxit = 1000
  
  full_fit <- #suppressMessages(
    emuFit_micro_penalized(X = X,
                           Y = Y,
                           B = NULL,
                           constraint_fn = rep(list(mean), 2), 
                           tolerance = 1e-7,
                           verbose = FALSE)
  
  B <- full_fit$B
  Y_aug <- full_fit$Y_augmented
  
  X_cup <- X_cup_from_X(X,J)

  j_ref <- 5
  null_fit <- fit_null(B = B,
                            Y = Y_aug,
                            X = X ,
                            X_cup = X_cup,
                            k_constr = k_constr,
                            j_constr = j_constr,
                            j_ref = j_ref,
                            constraint_fn = constraint_fn,
                            constraint_tol = 1e-5,
                            B_tol = 1e-4,
                            constraint_grad_fn =constraint_grad_fn,
                            verbose = FALSE,
                            trackB = FALSE) ## just track for one j

  
null_repar_fit <- fit_null_symmetric(Y = Y_aug,
                                X = X,
                                B = B,
                                j_constr = j_constr,
                                k_constr = k_constr,
                                j_ref = j_ref,
                                constraint_fn = constraint_fn[[1]],
                                constraint_grad_fn = constraint_grad_fn[[1]],
                                tolerance = 1e-4,
                                verbose = TRUE,
                                maxit = 1000)

  #min_mse_lag finds the lagrange multiplier that minimizes the squared norm
  #of the derivative of the lagrangian (ll + lambda*(g - B[k_constr,j_constr])
  #and returns the minimized squared norm -- lower means more accurate fit under 
  #null

  null_min_lag_norm <- min_mse_lag(X = X,
                                   Y = Y,
                                   B = null_fit$B,
                                   constraint_grad_fn = constraint_grad_fn[[1]],
                                   k_constr = k_constr,
                                   j_constr = j_constr,
                                   j_ref = j_ref)
  null_repar_min_lag_norm <- min_mse_lag(X = X,
                                   Y = Y,
                                   B = null_repar_fit$B,
                                   constraint_grad_fn = constraint_grad_fn[[1]],
                                   k_constr = k_constr,
                                   j_constr = j_constr,
                                   j_ref = j_ref)
  
  #sol'ns are at least close to equal
  expect_equal(null_repar_fit$B,null_fit$B,tolerance = 1e-2)
  
  #and to extent that they are not equal, repar fit is more accurate
  expect_true(null_min_lag_norm > null_repar_min_lag_norm)
              
        
})
