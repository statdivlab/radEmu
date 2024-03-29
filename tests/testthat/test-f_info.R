test_that("computing information matrix gives same result regardless of method", {

  Y <- matrix(1:4,nrow = 2)
  X <- cbind(1,c(0,1))
  penalized_fit <- emuFit_micro_penalized(X,
                                          Y,
                                          B = NULL,
                                          constraint_fn = function(x) x[2],
                                          maxit = 500,
                                          tolerance = 1e-5,
                                          verbose = FALSE)

  info1 <-
    f_info(Y,
         B_cup = B_cup_from_B(penalized_fit$B),
         B = penalized_fit$B,
         X = X,
         X_cup = X_cup_from_X(X,2),
         compute_together = TRUE)

  info2 <-
    f_info(Y,
           B_cup = B_cup_from_B(penalized_fit$B),
           B = penalized_fit$B,
           X = X,
           X_cup = X_cup_from_X(X,2),
           compute_together = FALSE)

  expect_equal(as.matrix(info1),as.matrix(info2))

})

test_that("Computed information is equal to numerical derivative with categorical predictor", {
  set.seed(59542234)
  n <- 2
  p <- 2
  X <- cbind(1,rep(c(0,1),each = n/2))
  J <- 2
  z <- rnorm(n) +3
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
  constraint_grad_fn <- function(x){mean(x)}

  ##### Arguments to fix:

  # constraint_grad_fn <- function(x){dpseudohuber_center_dx(x,0.1)
  constraint_grad_fn <- function(x){ rep(1/length(x), length(x))}
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
  I <- f_info(Y = Y_aug,
              B_cup = B_cup_from_B(B),
              B = B,
              X = X,
              X_cup = X_cup)

  ll <- function(b){
    tempB <- matrix(b, byrow = FALSE,ncol = 2)
    z <- update_z(Y= Y_aug,X = X, B = tempB)
    log_means <- do.call(cbind,
                         lapply(1:J,
                                function(j) z + X%*%tempB[,j,drop = FALSE]))
    return(sum(-Y_aug*log_means + exp(log_means)))
  }

  numo <- optim(rep(0,4),ll)

  optB <- matrix(numo$par, byrow = FALSE,ncol = 2)

  for(k in 1:2){
    optB[k,] <- optB[k,]- mean(optB[k,])
  }

  nh <- numDeriv::hessian(ll,c(B[,1],B[,2]))

  expect_true(max(abs(I - nh))<1e-5)

  # ll_repar <- function(b1){
  #   tempB <- cbind(b1,0)
  #   z <- update_z(Y= Y_aug,X = X, B = tempB)
  #   log_means <- do.call(cbind,
  #                        lapply(1:J,
  #                               function(j) z + X%*%tempB[,j,drop = FALSE]))
  #   return(sum(-Y_aug*log_means + exp(log_means)))
  # }
  #
  # numDeriv::hessian(ll_repar,c(B[,1]-B[,2]))


  #
  #
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
  #                j_ref = 1,
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

})
