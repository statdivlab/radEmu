test_that("Robust score statistic is invariant to reference taxon", {
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
                           tolerance = 1e-3,
                           verbose = FALSE)#)

  B <- full_fit$B
  Y_aug <- full_fit$Y_augmented



  X_cup <- X_cup_from_X(X,J)

  fit_thing <-
    score_test(B = B, #B (MPLE)
             Y = Y_aug, #Y (with augmentations)
             X = X, #design matrix
             X_cup = X_cup,
             k_constr = k_constr, #row index of B to constrain
             j_constr = j_constr,
             constraint_fn = constraint_fn, #constraint function
             constraint_grad_fn = constraint_grad_fn, #gradient of constraint fn
             rho_init = rho_init,
             tau = tau,
             kappa = kappa,
             j_ref = 5,
             B_tol = 1e-3,
             constraint_tol = constraint_tol,
             c1 = c1,
             maxit = maxit,
             inner_maxit = inner_maxit,
             verbose = FALSE)


  score_stats <-
    sapply(1:J,function(j_ref)
  get_score_stat(Y = Y_aug,
                 X_cup = X_cup,
                 X = X,
                 B = fit_thing$null_B,
                 k_constr = k_constr,
                 j_constr = j_constr,
                 constraint_grad_fn = constraint_grad_fn,
                 indexes_to_remove = (j_ref-1)*p + 1:p ,
                 j_ref = j_ref,
                 J = J,
                 n = n,
                 p = p,
                 I = NULL,
                 Dy = NULL)
    )

  expect_true(sd(score_stats)<1e-5)


})
#


test_that("Sandwich score statistic is asymptotically chi square distributed
under null when Poisson assumption is met", {

  skip("Skipping test that requires 1000 simulations be run")
  set.seed(595434)
  n <- 10
  p <- 2
  X <- cbind(1,rep(c(0,1),each = n/2))
  J <- 25
  z <- rnorm(n) +8
  b0 <- rnorm(J)
  b1 <- seq(1,10,length.out = J)
  # b1 <- c(-3,-2,0,5,0)
  b1 <- b1 - mean(b1)
  # b0 <- b0 - mean(b0)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)

  k_constr <- 2
  j_constr <- 13
  p <- 2

  constraint_fn <- function(x){ pseudohuber_center(x,0.1)}
  # constraint_fn <- function(x){mean(x)}

  ##### Arguments to fix:

  constraint_grad_fn <- function(x){dpseudohuber_center_dx(x,0.1)}
  # constraint_grad_fn <- function(x){rep(1/length(x),length(x))}
  rho_init = 1
  tau = 2
  kappa = 0.8
  obj_tol = 100
  B_tol <- 1e-3
  constraint_tol = 1e-5
  init_tol = 1e6
  c1 = 1e-4
  maxit = 1000
  inner_maxit = 25

  nsim <- 1000
  score_stats <- numeric(nsim)
  for(sim in 1:nsim){
    Y[] <- 0
    for(i in 1:n){
      while(sum(Y[i,])==0){
        for(j in 1:J){
          temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
          # Y[i,j] <- rpois(1, lambda = temp_mean)
          Y[i,j] <- rnbinom(1,mu = temp_mean, size = 3)*rbinom(1,1,0.6)
        }
      }
    }

    full_fit <- #suppressMessages(
      emuFit_micro_penalized(X = X,
                             Y = Y,
                             B = NULL,
                             constraint_fn = mean,
                             tolerance = 1e-3,
                             verbose = FALSE)#)

    B <- full_fit$B
    Y_aug <- full_fit$Y_augmented



    X_cup <- X_cup_from_X(X,J)

    fit_thing <-
      score_test(B = B, #B (MPLE)
                 Y = Y_aug, #Y (with augmentations)
                 X = X, #design matrix
                 X_cup = X_cup,
                 k_constr = k_constr, #row index of B to constrain
                 j_constr = j_constr,
                 constraint_fn = constraint_fn, #constraint function
                 constraint_grad_fn = constraint_grad_fn, #gradient of constraint fn
                 rho_init = rho_init,
                 tau = tau,
                 kappa = kappa,
                 j_ref = 5,
                 B_tol = 0.001,
                 constraint_tol = constraint_tol,
                 c1 = c1,
                 maxit = maxit,
                 inner_maxit = inner_maxit,
                 verbose = FALSE)

    j_ref <- 5

    score_stats[sim] <- get_score_stat(Y = Y_aug,
                                       X_cup = X_cup,
                                       X = X,
                                       B = fit_thing$null_B,
                                       k_constr = k_constr,
                                       j_constr = j_constr,
                                       constraint_grad_fn = constraint_grad_fn,
                                       indexes_to_remove = (j_ref-1)*p + 1:p ,
                                       j_ref = 5,
                                       J = J,
                                       n = n,
                                       p = p,
                                       check_influence = FALSE,
                                       I = NULL,
                                       Dy = NULL,
                                       model_based = FALSE)

    x <- seq(-10,5,0.01)
    # hist(log(score_stats[1:sim]),breaks = 20,freq = FALSE)
    # lines(x,exp(x)*dchisq(,1),lty = 2, col = "red")
    e_cdf <- sapply(x,function(y) mean(score_stats[1:sim]<=exp(y)))
    plot(x,e_cdf,type = "l")
    lines(x, e_cdf + sqrt(e_cdf*(1 - e_cdf)/sim)*2,lty = 2)
    lines(x, e_cdf - sqrt(e_cdf*(1 - e_cdf)/sim)*2,lty = 2)
    lines(x,pchisq(exp(x),1),lty = 3, col = "red")
    text(0,0.2,sim)
    text(3,0.2,round(mean(score_stats[1:sim]),2))
  }



})

test_that("model-based score statistic is invariant to reference taxon", {
  set.seed(59542234)
  n <- 10
  p <- 2
  X <- cbind(1,rep(c(0,1),each = n/2))
  J <- 5
  z <- rnorm(n) +8
  b0 <- rnorm(J)
  b1 <- seq(1,3,length.out = J)
  # b1 <- c(-3,-2,0,5,0)
  b1 <- b1 - mean(b1)
  # b0 <- b0 - mean(b0)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)

  k_constr <- 2
  j_constr <- 3
  p <- 2

  # constraint_fn <- function(x){ pseudohuber_center(x,0.1)}
  constraint_fn <- function(x){mean(x)}

  ##### Arguments to fix:

  # constraint_grad_fn <- function(x){dpseudohuber_center_dx(x,0.1)
  constraint_grad_fn <- function(x){rep(1/length(x),length(x))}
  rho_init = 1
  tau = 2
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
                             tolerance = 1e-3,
                             verbose = FALSE)#)

    B <- full_fit$B
    Y_aug <- full_fit$Y_augmented



    X_cup <- X_cup_from_X(X,J)

    fit_thing <-
      score_test(B = B, #B (MPLE)
                 Y = Y_aug, #Y (with augmentations)
                 X = X, #design matrix
                 X_cup = X_cup,
                 k_constr = k_constr, #row index of B to constrain
                 j_constr = j_constr,
                 constraint_fn = constraint_fn, #constraint function
                 constraint_grad_fn = constraint_grad_fn, #gradient of constraint fn
                 rho_init = rho_init,
                 tau = tau,
                 kappa = kappa,
                 j_ref = 5,
                 B_tol = 0.0001,
                 constraint_tol = constraint_tol,
                 c1 = c1,
                 maxit = maxit,
                 inner_maxit = inner_maxit,
                 verbose = FALSE)


    score_stats <- sapply(1:J,function(jref){
      get_score_stat(Y = Y_aug,
                                       X_cup = X_cup,
                                       X = X,
                                       B = fit_thing$null_B,
                                       k_constr = k_constr,
                                       j_constr = j_constr,
                                       constraint_grad_fn = constraint_grad_fn,
                                       indexes_to_remove = (jref-1)*p + 1:p ,
                                       j_ref = jref,
                                       J = J,
                                       n = n,
                                       p = p,
                                       I = NULL,
                                       Dy = NULL)
    })

   expect_true(sd(score_stats)<1e-8)

})

#
# test_that("multiplication works", {
#   set.seed(59542234)
#   n <- 10
#   p <- 2
#   X <- cbind(1,rep(c(0,1),each = n/2))
#   J <- 2
#   z <- rnorm(n) +3
#   b0 <- rnorm(J)
#   b1 <- seq(1,10,length.out = J)
#   b1 <- b1 - mean(b1)
#   b0 <- b0 - mean(b0)
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = J, nrow = n)
#
#   k_constr <- 2
#   j_constr <- 1
#   p <- 2
#
#   # constraint_fn <- function(x){ pseudohuber_center(x,0.1)}
#   constraint_grad_fn <- function(x){x[2]}
#
#   ##### Arguments to fix:
#
#   # constraint_grad_fn <- function(x){dpseudohuber_center_dx(x,0.1)
#   constraint_grad_fn <- function(x){ 0:1}
#   rho_init = 1
#   tau = 1.2
#   kappa = 0.8
#   obj_tol = 100
#   score_tol <- 1e-3
#   constraint_tol = 1e-5
#   init_tol = 1e6
#   c1 = 1e-4
#   maxit = 1000
#   inner_maxit = 25
#
#   Y[] <- 0
#   for(i in 1:n){
#     while(sum(Y[i,])==0){
#       for(j in 1:J){
#         temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#         Y[i,j] <- rpois(1, lambda = temp_mean)
#         # Y[i,j] <- rnbinom(1,mu = temp_mean, size = 3)*rbinom(1,1,0.6)
#       }
#     }
#   }
#
#   full_fit <- #suppressMessages(
#     emuFit_micro_penalized(X = X,
#                            Y = Y,
#                            B = NULL,
#                            constraint_fn = mean,
#                            tolerance = 1e-3,
#                            verbose = TRUE)#)
#
#   B <- full_fit$B
#   Y_aug <- full_fit$Y_augmented
#
#
#
#   X_cup <- X_cup_from_X(X,J)
#
#   fit_thing <-
#     score_test(B = B, #B (MPLE)
#                Y = Y_aug, #Y (with augmentations)
#                X = X, #design matrix
#                X_cup = X_cup,
#                k_constr = k_constr, #row index of B to constrain
#                j_constr = j_constr,
#                constraint_fn = constraint_fn, #constraint function
#                constraint_grad_fn = constraint_grad_fn, #gradient of constraint fn
#                rho_init = rho_init,
#                tau = tau,
#                kappa = kappa,
#                j_ref = J,
#                score_tol = 5,
#                constraint_tol = constraint_tol,
#                c1 = c1,
#                maxit = maxit,
#                inner_maxit = inner_maxit,
#                verbose = TRUE)
#
#
#   score_stats <-
#     sapply(1:J,function(j_ref)
#       get_score_stat(Y = Y_aug,
#                      X_cup = X_cup,
#                      X = X,
#                      B = fit_thing$null_B,
#                      k_constr = 2,
#                      j_constr = 1,
#                      constraint_grad_fn = constraint_grad_fn,
#                      indexes_to_remove = (j_ref-1)*p + 1:p ,
#                      j_ref = j_ref,
#                      J = J,
#                      n = n,
#                      p = p,
#                      check_influence = FALSE,
#                      I = NULL,
#                      Dy = NULL)
#     )
#
#   plot(score_stats)
#   expect_true(max(score_stats) - min(score_stats)<1e-3)
#
# })
#
#
