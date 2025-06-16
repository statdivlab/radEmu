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
  
  # constraint_fn <- function(x){ pseudohuber_median(x,0.1)}
  constraint_fn <- rep(list(function(x){mean(x)}), 2)
  
  ##### Arguments to fix:
  
  # constraint_grad_fn <- function(x){dpseudohuber_median_dx(x,0.1)
  constraint_grad_fn <- rep(list(function(x){ rep(1/length(x), length(x))}), 2)
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
  
  full_fit <- #suppressMessages(
    emuFit_micro_penalized(X = X,
                           Y = Y,
                           B = NULL,
                           constraint_fn = rep(list(mean), 2), 
                           tolerance = 1e-7,
                           verbose = FALSE)
  
  B <- full_fit$B
  Y_aug <- full_fit$Y_augmented
  
  X_cup <- X_cup_from_X_fast(X,J)
  
  fits <- vector(5, mode = "list")
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
                              verbose = FALSE,
                              trackB = (j_ref == 4)) ## just track for one j
  }
  
  Bs <- lapply(fits,function(x){
    temp_B <- x$B
    for(k in 1:2){
      temp_B[k,] <- temp_B[k,] - constraint_fn[[k]](temp_B[k,])
    }
    return(B_cup_from_B(temp_B))
  })
  
  ## check results are same regardless of column index of convenience constraint
  expect_equal(Bs[[1]], Bs[[2]],tolerance = 1e-2)
  expect_equal(Bs[[2]], Bs[[3]],tolerance = 1e-2)
  expect_equal(Bs[[3]], Bs[[4]],tolerance = 1e-2)
  expect_equal(Bs[[4]], Bs[[5]],tolerance = 1e-2)
  expect_equal(mean(fits[[1]]$B[2,]), fits[[1]]$B[2,1],tolerance = 1e-4)
  
  # Null fit satisfies null constraints. 
  expect_equal(Bs[[1]][k_constr,j_constr], constraint_fn[[k_constr]](Bs[[1]][k_constr]))
  expect_equal(Bs[[2]][k_constr,j_constr], constraint_fn[[k_constr]](Bs[[2]][k_constr]))
  expect_equal(Bs[[4]][k_constr,j_constr], constraint_fn[[k_constr]](Bs[[4]][k_constr]))
  
  # If trackB = TRUE we get B at each iteration back
  expect_true(max(fits[[4]]$Bs$iter) > 20)
  
})

