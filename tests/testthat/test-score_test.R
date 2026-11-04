test_that("We get same score test results regardless of whether we provide inverse of info *and*
we do *not* get same results if we use incorrect info", {
  
  set.seed(343234)
  J <- 10
  n <- 40
  X <- cbind(1,rep(c(0,1),each = n/2))
  b0 <- rnorm(J)
  b1 <- 1:J
  b1 <- b1 - mean(b1)
  b1[5] <- pseudohuber_median(b1[-5],0.1)
  b0 <- b0 - mean(b0)
  b <- rbind(b0, b1)
  Y <- radEmu:::simulate_data(n = n,
                              J = J,
                              X = X,
                              b0 = b0,
                              b1 = b1,
                              distn = "ZINB",
                              zinb_size = 3,
                              zinb_zero_prop = 0.6,
                              mean_z = 4)

  k_constr <- 2
  j_constr <- 5
  p <- 2



  constraint_fn <- function(x){ pseudohuber_median(x,0.1)}

  ##### Arguments to fix:

    full_fit <- emuFit(X = X,
                       Y = Y,
                       B = NULL,
                       tolerance = 0.01,
                       verbose = FALSE,
                       run_score_tests = FALSE,
                       match_row_names = FALSE)


    B <- full_fit$B
    Y_aug <- full_fit$Y_augmented

    j_ref <- which.max(colSums(Y>0))

    indexes_to_remove <- (j_ref - 1)*p + 1:p


    X_cup <- X_cup_from_X(X,J)

    score_test_as_is <-
      score_test(B = B,
                 Y = Y_aug,
                 X = X,
                 X_cup = X_cup,
                 constraint_fn = rep(list(function(x) pseudohuber_median(x,0.1)), 2), 
                 constraint_grad_fn = rep(list(function(x) dpseudohuber_median_dx(x,0.1)), 2), 
                 j_ref = j_ref,
                 B_tol = 0.01,
                 constraint_tol= 0.001,
                 verbose = FALSE,
                 k_constr = 2,
                 j_constr = 1)

    I <- f_info(Y = Y_aug,
                B_cup = B_cup_from_B(score_test_as_is$null_B),
                B = score_test_as_is$null_B,
                X = X,
                X_cup = X_cup)

    I_inv <- Matrix::solve(I[-indexes_to_remove,-indexes_to_remove])

    score_stat_with_I_inv <-
      get_score_stat(Y = Y_aug,
                     X = X,
                     X_cup = X_cup,
                     B =    score_test_as_is$null_B,
                     k_constr = 2,
                     j_constr = 1,
                     constraint_grad_fn = rep(list(function(x) dpseudohuber_median_dx(x,0.1)), 2), 
                     indexes_to_remove = indexes_to_remove,
                     j_ref = j_ref,
                     J = J,
                     n = n,
                     p = p,
                     I_inv = I_inv)$score_stat

    score_stat_with_other_mat <-
      get_score_stat(Y = Y_aug,
                     X = X,
                     X_cup = X_cup,
                     B =    score_test_as_is$null_B,
                     k_constr = 2,
                     j_constr = 1,
                     constraint_grad_fn = rep(list(function(x) dpseudohuber_median_dx(x,0.1)), 2), 
                     indexes_to_remove = indexes_to_remove,
                     j_ref = j_ref,
                     J = J,
                     n = n,
                     p = p,
                     I_inv = diag(rep(1,18)))$score_stat

    expect_equal(score_stat_with_I_inv,score_test_as_is$score_stat)
    expect_true(score_stat_with_other_mat !=score_test_as_is$score_stat)

})

# Small data set for the "fit_null errors" tests below. The values don't matter,
# since those tests only check how a failing `fit_null()` is handled.
make_fit_null_error_data <- function() {
  set.seed(1)
  n <- 10
  J <- 5
  X <- cbind(1, rep(0:1, each = n / 2))
  Y <- matrix(rpois(n * J, 5), n, J)
  list(X = X, Y = Y, J = J)
}

# `score_test()` wraps `fit_null()` in a non-silent `try()`, which prints each
# caught error with `cat()` to `getOption("try.outFile")` (stderr by default).
# That isn't a message or warning, so `suppressMessages()` and
# `suppressWarnings()` can't silence it. Send that output to the null device
# while `expr` runs. Warnings still propagate, so `expect_warning()` works.
without_try_output <- function(expr) {
  old <- options(try.outFile = nullfile())
  on.exit(options(old))
  expr
}

test_that("if fit_null errors on every attempt, score_test warns and returns NULL", {
  local_mocked_bindings(
    fit_null = function(...) stop("injected fit_null failure")
  )
  
  dat <- make_fit_null_error_data()
  X <- dat$X
  Y <- dat$Y
  J <- dat$J

  # `fit_null()` is mocked, but the full (unconstrained) fit doesn't use it.
  fit <- suppressMessages(
    emuFit(Y = Y, X = X, run_score_tests = FALSE, verbose = FALSE)
  )

  # `constraint_fn` has no "constraint_type" attribute, so `score_test()` treats
  # it as "other" and goes straight to the `fit_null()` retry loop. It needs one
  # entry per row of B (p = 2) because it is indexed by `k_constr`.
  expect_warning(
    res <- without_try_output(suppressMessages(score_test(
      B = fit$B,
      Y = fit$Y_augmented,
      X = X,
      X_cup = X_cup_from_X_fast(X, J),
      k_constr = 2,
      j_constr = 1,
      j_ref = J,
      ntries = 2,
      constraint_fn = rep(list(mean), 2),
      constraint_grad_fn = rep(
        list(function(x) rep(1 / length(x), length(x))),
        2
      )
    ))),
    "failed with an error \\(injected fit_null failure\\)"
  )
  
  expect_null(res)
})

test_that("if fit_null errors on every attempt, emuFit reports NA and marks the tests as failed", {
  local_mocked_bindings(
    fit_null = function(...) stop("injected fit_null failure")
  )
  
  dat <- make_fit_null_error_data()
  X <- dat$X
  Y <- dat$Y

  # The mocked `fit_null()` is only reached when constraint_type is "other",
  # which is what `null_fit_alg = "augmented_lagrangian"` sets. With J = 5,
  # emuFit already forces augmented Lagrangian (it does so whenever J < 20), so
  # the argument is technically redundant. It's kept explicit so the test still
  # works if J is raised: for J >= 20, this X (intercept + one binary covariate)
  # would make emuFit pick the "discrete" algorithm, which calls
  # `fit_null_discrete()` and might never reach `fit_null()`.
  #
  # `ntries = 2` (default 4) just makes the test faster; every attempt fails.
  full <- without_try_output(suppressWarnings(suppressMessages(emuFit(
    Y = Y,
    X = X,
    test_kj = data.frame(k = 2, j = 1:2),
    null_fit_alg = "augmented_lagrangian",
    control = list(ntries = 2),
    verbose = FALSE
  ))))
  
  expect_true(all(is.na(full$coef$score_stat[1:2])))
  expect_true(all(is.na(full$coef$pval[1:2])))
  expect_equal(
    full$score_test_hyperparams$converged,
    rep("failed with error", 2)
  )
  expect_equal(full$null_estimation_unconverged$j, 1:2)
})

# test_that("simple negative binomial example fits reasonably fast and returns reasonable output", {
#   
#   set.seed(343234)
#   X <- cbind(1,rep(c(0,1),each = 20))
#   J <- 10
#   n <- 40
#   z <- rnorm(40) +8
#   b0 <- rnorm(10)
#   b1 <- 1:10
#   b1 <- b1 - mean(b1)
#   b1[5] <- pseudohuber_median(b1[-5],0.1)
#   b0 <- b0 - mean(b0)
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 10, nrow = 40)
#
#   k_constr <- 2
#   j_constr <- 5
#   p <- 2
#
#
#
#   constraint_fn <- function(x){ pseudohuber_median(x,0.1)}
#
#   ##### Arguments to fix:
#
#   constraint_grad_fn <- function(x){dpseudohuber_median_dx(x,0.1)}
#
#   constraint_hess_fn <- function(x,ind_1,ind_2){hess_pseudohuber_median(x,0.1,ind_1,ind_2)}
#   rho_init = 1
#   tau = 10
#   kappa = 0.8
#   obj_tol = 1e-6
#   constraint_tol = 1e-5
#   init_tol = 100
#   c1 = 1e-4
#   maxit = 1000
#   inner_maxit = 25
#   nsim <- 10
#   score_stats <- vector(nsim,mode = "list")
#   for(sim in 1:nsim){
#     message("Simulation ", sim, " of ",nsim,".")
#   for(i in 1:40){
#     for(j in 1:10){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       # Y[i,j] <- rpois(1, lambda = temp_mean)
#       Y[i,j] <- rnbinom(1,mu = temp_mean, size = 3)*rbinom(1,1,0.6)
#     }
#   }
#
#   full_fit <- suppressMessages(emuFit_micro_penalized(X = X,
#                                                       Y = Y,
#                                                       B = NULL,
#                                                       constraint_fn = mean,
#                                                       verbose = TRUE))
#
#   B <- full_fit$B
#   Y_aug <- full_fit$Y_augmented

#   X_cup = X_cup_from_X(X,J)
#   scores <- vector(n,mode = "list")
#
#   #compute score contributions of observations i = 1 through n
#   for(i in 1:n){
#     X_cup_i <- X_cup[(i - 1)*J + 1:J,]
#     scores[[i]] <- as.matrix(dpll_dB_cup(X[i,,drop = FALSE],Y[i,,drop = FALSE],B))
#   }
#   Dy <- Reduce("+",lapply(scores,function(x) tcrossprod(x)))
#   Y <- Y_aug
#
#   score_results <-
#     lapply(exp(seq(log(1e4),log(1e-8),length.out = 13)), function(ot){
#       print(ot);
#     score_test(B = B, #B (MPLE)
#                Y = Y, #Y (with augmentations)
#                X = X, #design matrix
#                k_constr = k_constr, #row index of B to constrain
#                j_constr = j_constr,
#                constraint_fn = constraint_fn, #constraint function
#                constraint_grad_fn = constraint_grad_fn, #gradient of constraint fn
#                # constraint_hess_fn = constraint_hess_fn, #hessian of constraint fn
#                rho_init = rho_init,
#                tau = tau,
#                kappa = kappa,
#                obj_tol = ot,
#                constraint_tol = constraint_tol,
#                init_tol = init_tol,
#                c1 = c1,
#                maxit = maxit,
#                inner_maxit = inner_maxit,
#                Dy = Dy,
#                verbose = FALSE)})
#
#   score_stats[[sim]] <- sapply(1:length(score_results),function(i) score_results[[i]]$score_stat)
#   }

#   do.call(rbind,lapply(1:nsim,function(x)
#     data.frame(log_tol = seq(log(1e4),log(1e-8),length.out = 13)/log(10),
#                score_stat = score_stats[[x]],
#                simulation = x))) %>%
#     group_by(simulation) %>%
#     mutate(score_stat = score_stat - score_stat[which.min(log_tol)]) %>%
#
#     ungroup %>%
#     ggplot() +
#     geom_point(aes(x = log_tol,y = score_stat,color = as.factor(simulation))) +
#     geom_line(aes(x = log_tol,y = score_stat,color = as.factor(simulation),
#                   group = as.factor(simulation))) +
#     theme_bw()
#
#
#
# })
#
#
# test_that("simple negative binomial example fits reasonably fast and returns reasonable output", {
#   
#   set.seed(343234)
#   n <- 250
#   X <- cbind(1,rep(c(0,1),each = n/2))
#   J <- 10
#   z <- rnorm(n) +8
#   b0 <- rnorm(10)
#   b1 <- 1:10
#   b1 <- b1 - mean(b1)
#   b1[5] <- pseudohuber_median(b1[-5],0.1)
#   b0 <- b0 - mean(b0)
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 10, nrow = n)
#
#   k_constr <- 2
#   j_constr <- 5
#   p <- 2
#
#
#
#   constraint_fn <- function(x){ pseudohuber_median(x,0.1)}
#
#   ##### Arguments to fix:
#
#   constraint_grad_fn <- function(x){dpseudohuber_median_dx(x,0.1)}
#
#   constraint_hess_fn <- function(x,ind_1,ind_2){hess_pseudohuber_median(x,0.1,ind_1,ind_2)}
#   rho_init = 1
#   tau = 10
#   kappa = 0.8
#   obj_tol = 1e-6
#   constraint_tol = 1e-5
#   init_tol = 1
#   c1 = 1e-4
#   maxit = 1000
#   inner_maxit = 25
#   nsim <- 1000
#   score_stats <- numeric(nsim)
#   for(sim in 1:nsim){
#     message("Simulation ", sim, " of ",nsim,".")
#     Y[] <- 0
#     for(i in 1:n){
#       while(sum(Y[i,])==0){
#       for(j in 1:10){
#         temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#         # Y[i,j] <- rpois(1, lambda = temp_mean)
#         Y[i,j] <- rnbinom(1,mu = temp_mean, size = 3)*rbinom(1,1,0.6)
#       }
#         }
#     }
#
#     full_fit <- suppressMessages(
#       emuFit_micro_penalized(X = X,
#                              Y = Y,
#                              B = NULL,
#                              constraint_fn = mean,
#                              verbose = TRUE))
#
#     B <- full_fit$B
#     Y_aug <- full_fit$Y_augmented
#
#     X_cup = X_cup_from_X(X,J)
#     scores <- vector(n,mode = "list")
#
#     #compute score contributions of observations i = 1 through n
#     for(i in 1:n){
#       X_cup_i <- X_cup[(i - 1)*J + 1:J,]
#       scores[[i]] <- as.matrix(dpll_dB_cup(X[i,,drop = FALSE],Y[i,,drop = FALSE],B))
#     }
#     Dy <- Reduce("+",lapply(scores,function(x) tcrossprod(x)))
#     Y <- Y_aug
#
#     score_results <-
#         score_test(B = B, #B (MPLE)
#                    Y = Y, #Y (with augmentations)
#                    X = X, #design matrix
#                    k_constr = k_constr, #row index of B to constrain
#                    j_constr = j_constr,
#                    constraint_fn = constraint_fn, #constraint function
#                    constraint_grad_fn = constraint_grad_fn, #gradient of constraint fn
#                    constraint_hess_fn = constraint_hess_fn, #hessian of constraint fn
#                    rho_init = rho_init,
#                    tau = tau,
#                    kappa = kappa,
#                    obj_tol = obj_tol,
#                    constraint_tol = constraint_tol,
#                    init_tol = init_tol,
#                    c1 = c1,
#                    maxit = maxit,
#                    inner_maxit = inner_maxit,
#                    Dy = Dy,
#                    verbose = TRUE)
#
#     score_stats[sim] <- score_results$score_stat
#
#     if(sim%%10==0){
#       type1 <- mean(pchisq(score_stats[1:sim],1,lower.tail = FALSE)<=0.05)
#       print(type1)
#       x <- exp(seq(-10,3,0.1))
#       e_cdf <- sapply(x,function(k) mean(score_stats[1:sim]<=k))
#       plot(log(x),e_cdf,type = "s",
#            xlab = "score statistic",ylab = "cumulative probability",
#            main = paste("Summary of first",sim,"simulations for n =",n),
#            sub = paste(signif(type1,3)))
#       lines(log(x),e_cdf + (2/sqrt(sim))*sqrt(e_cdf*(1 - e_cdf)),lty = 3,type = "s")
#       lines(log(x),e_cdf - (2/sqrt(sim))*sqrt(e_cdf*(1 - e_cdf)),lty = 3,type = "s")
#       # hist(score_stats[1:sim],breaks = 50,freq = FALSE)
#       lines(log(x),pchisq(x,1),lty = 2, col = "red")
#     }
#   }
#
#
#   hist(score_stats[1:75],breaks = 20,freq = FALSE)
#   lines(x = seq(0,10,0.001),dchisq(seq(0,10,0.001),1),lty= 2)
#
#   # pchisq(score_stats[1:75],1,lower.tail = FALSE) %>%
#
#
#
# })


#
# test_that("simple negative binomial example fits reasonably fast and returns reasonable output", {
#   
#   set.seed(343234)
#   n <- 10
#   X <- cbind(1,rep(c(0,1),each = n/2))
#   J <- 10
#   z <- rnorm(n) +8
#   b0 <- rnorm(10)
#   b1 <- 1:10
#   b1 <- b1 - mean(b1)
#   b1[5] <- pseudohuber_median(b1[-5],0.1)
#   b0 <- b0 - mean(b0)
#   b <- rbind(b0,b1)
#   Y <- matrix(0,ncol = 10, nrow = n)
#
#   k_constr <- 2
#   j_constr <- 5
#   p <- 2
#
#
#
#   constraint_fn <- function(x){ pseudohuber_median(x,0.1)}
#
#   ##### Arguments to fix:
#
#   constraint_grad_fn <- function(x){dpseudohuber_median_dx(x,0.1)}
#
#   constraint_hess_fn <- function(x,ind_1,ind_2){hess_pseudohuber_median(x,0.1,ind_1,ind_2)}
#   rho_init = 1
#   tau = 10
#   kappa = 0.8
#   obj_tol = 1e-5
#   constraint_tol = 1e-5
#   init_tol = 1
#   c1 = 1e-4
#   maxit = 1000
#   inner_maxit = 25
#   nsim <- 1000
#   score_stats <- numeric(nsim)
#   for(sim in 1:nsim){
#     message("Simulation ", sim, " of ",nsim,".")
#     Y[] <- 0
#     for(i in 1:n){
#       while(sum(Y[i,])==0){
#       for(j in 1:10){
#         temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#         # Y[i,j] <- rpois(1, lambda = temp_mean)
#         Y[i,j] <- rnbinom(1,mu = temp_mean, size = 3)*rbinom(1,1,0.6)
#       }
#         }
#     }
#
#     full_fit <- suppressMessages(
#       emuFit_micro_penalized(X = X,
#                              Y = Y,
#                              B = NULL,
#                              constraint_fn = mean,
#                              verbose = TRUE))
#
#     B <- full_fit$B
#     Y_aug <- full_fit$Y_augmented
#
#
#     X_cup = X_cup_from_X(X,J)
#     scores <- vector(n,mode = "list")
#
#     #compute score contributions of observations i = 1 through n
#     for(i in 1:n){
#       X_cup_i <- X_cup[(i - 1)*J + 1:J,]
#       scores[[i]] <- as.matrix(dpll_dB_cup(X[i,,drop = FALSE],Y_aug[i,,drop = FALSE],B))
#     }
#     Dy <- Reduce("+",lapply(scores,function(x) tcrossprod(x)))
#     old_Y <- Y
#     Y <- Y_aug
#
#     score_results <-
#       score_test(B = B, #B (MPLE)
#                  Y = Y, #Y (with augmentations)
#                  X = X, #design matrix
#                  k_constr = k_constr, #row index of B to constrain
#                  j_constr = j_constr,
#                  constraint_fn = constraint_fn, #constraint function
#                  constraint_grad_fn = constraint_grad_fn, #gradient of constraint fn
#                  constraint_hess_fn = constraint_hess_fn, #hessian of constraint fn
#                  rho_init = rho_init,
#                  tau = tau,
#                  kappa = kappa,
#                  obj_tol = obj_tol,
#                  constraint_tol = constraint_tol,
#                  init_tol = init_tol,
#                  c1 = c1,
#                  maxit = maxit,
#                  inner_maxit = inner_maxit,
#                  Dy = Dy,
#                  verbose = TRUE)
#
#     score_stats[sim] <- score_results$score_stat
#
#
#     if(sim%%10==0){
#       type1 <- mean(pchisq(score_stats[1:sim],1,lower.tail = FALSE)<=0.05)
#       print(type1)
#       x <- exp(seq(-10,3,0.1))
#       e_cdf <- sapply(x,function(k) mean(score_stats[1:sim]<=k))
#       plot(log(x),e_cdf,type = "s",
#            xlab = "score statistic",ylab = "cumulative probability",
#            main = paste("Summary of first",sim,"simulations for n =",n),
#            sub = paste(signif(type1,3)))
#       lines(log(x),e_cdf + (2/sqrt(sim))*sqrt(e_cdf*(1 - e_cdf)),lty = 3,type = "s")
#       lines(log(x),e_cdf - (2/sqrt(sim))*sqrt(e_cdf*(1 - e_cdf)),lty = 3,type = "s")
#       # hist(score_stats[1:sim],breaks = 50,freq = FALSE)
#       lines(log(x),pchisq(x,1),lty = 2, col = "red")
#     }
#   }
#
#
#   hist(score_stats[1:75],breaks = 20,freq = FALSE)
#   lines(x = seq(0,10,0.001),dchisq(seq(0,10,0.001),1),lty= 2)
#
#   # pchisq(score_stats[1:75],1,lower.tail = FALSE) %>%
#
#
#
# })
#
# test_that("simple negative binomial example fits reasonably fast and returns reasonable output", {
#   
#   set.seed(343234)
#   n <- 50
#   X <- cbind(1,rep(c(0,1),each = n/2))
#   J <- 10
#   z <- rnorm(n) +8
#   b0 <- rnorm(10)
#   b1 <- 1:10
#   b1 <- b1 - mean(b1)
#   b1[5] <- pseudohuber_median(b1[-5],0.1)
#   b0 <- b0 - mean(b0)
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 10, nrow = n)
#
#   k_constr <- 2
#   j_constr <- 5
#   p <- 2
#
#
#
#   constraint_fn <- function(x){ pseudohuber_median(x,0.1)}
#
#   ##### Arguments to fix:
#
#   constraint_grad_fn <- function(x){dpseudohuber_median_dx(x,0.1)}
#
#   constraint_hess_fn <- function(x,ind_1,ind_2){hess_pseudohuber_median(x,0.1,ind_1,ind_2)}
#   rho_init = 1
#   tau = 10
#   kappa = 0.8
#   obj_tol = 1e-6
#   constraint_tol = 1e-5
#   init_tol = 1
#   c1 = 1e-4
#   maxit = 1000
#   inner_maxit = 25
#   nsim <- 1000
#   score_stats <- numeric(nsim)
#   for(sim in 1:nsim){
#     message("Simulation ", sim, " of ",nsim,".")
#     for(i in 1:n){
#       for(j in 1:10){
#         temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#         # Y[i,j] <- rpois(1, lambda = temp_mean)
#         Y[i,j] <- rnbinom(1,mu = temp_mean, size = 3)*rbinom(1,1,0.6)
#       }
#     }
#
#     full_fit <- suppressMessages(
#       emuFit_micro_penalized(X = X,
#                              Y = Y,
#                              B = NULL,
#                              constraint_fn = mean,
#                              verbose = TRUE))
#
#     B <- full_fit$B
#     Y_aug <- full_fit$Y_augmented
#
#     X_cup = X_cup_from_X(X,J)
#     scores <- vector(n,mode = "list")
#
#     #compute score contributions of observations i = 1 through n
#     for(i in 1:n){
#       X_cup_i <- X_cup[(i - 1)*J + 1:J,]
#       scores[[i]] <- as.matrix(dpll_dB_cup(X[i,,drop = FALSE],Y[i,,drop = FALSE],B))
#     }
#     Dy <- Reduce("+",lapply(scores,function(x) tcrossprod(x)))
#     Y <- Y_aug
#
#     score_results <-
#       score_test(B = B, #B (MPLE)
#                  Y = Y, #Y (with augmentations)
#                  X = X, #design matrix
#                  k_constr = k_constr, #row index of B to constrain
#                  j_constr = j_constr,
#                  constraint_fn = constraint_fn, #constraint function
#                  constraint_grad_fn = constraint_grad_fn, #gradient of constraint fn
#                  # constraint_hess_fn = constraint_hess_fn, #hessian of constraint fn
#                  rho_init = rho_init,
#                  tau = tau,
#                  kappa = kappa,
#                  obj_tol = obj_tol,
#                  constraint_tol = constraint_tol,
#                  init_tol = init_tol,
#                  c1 = c1,
#                  maxit = maxit,
#                  inner_maxit = inner_maxit,
#                  Dy = Dy,
#                  verbose = TRUE)
#
#     score_stats[sim] <- score_results$score_stat
#
#     if(sim%%10==0){
#       type1 <- mean(pchisq(score_stats[1:sim],1,lower.tail = FALSE)<=0.05)
#       print(type1)
#       x <- exp(seq(-10,3,0.1))
#       e_cdf <- sapply(x,function(k) mean(score_stats[1:sim]<=k))
#       plot(log(x),e_cdf,type = "s",
#            xlab = "score statistic",ylab = "cumulative probability",
#            main = paste("Summary of first",sim,"simulations for n =",n),
#            sub = paste(signif(type1,3)))
#       lines(log(x),e_cdf + (2/sqrt(sim))*sqrt(e_cdf*(1 - e_cdf)),lty = 3,type = "s")
#       lines(log(x),e_cdf - (2/sqrt(sim))*sqrt(e_cdf*(1 - e_cdf)),lty = 3,type = "s")
#       # hist(score_stats[1:sim],breaks = 50,freq = FALSE)
#       lines(log(x),pchisq(x,1),lty = 2, col = "red")
#     }
#   }
#
#
#   hist(score_stats[1:75],breaks = 20,freq = FALSE)
#   lines(x = seq(0,10,0.001),dchisq(seq(0,10,0.001),1),lty= 2)
#
#   # pchisq(score_stats[1:75],1,lower.tail = FALSE) %>%
#
#
#
# })
