test_that("we get same null fit with different j_ref", {
  set.seed(59542234)
  n <- 10
  J <- 5
  X <- cbind(1, rep(c(0, 1), each = n / 2))
  b0 <- rnorm(J)
  b1 <- seq(1, 10, length.out = J)
  b1 <- b1 - mean(b1)
  b0 <- b0 - mean(b0)
  Y <- radEmu:::simulate_data(
    n = n,
    J = J,
    X = X,
    b0 = b0,
    b1 = b1,
    distn = "Poisson",
    mean_z = 8
  )

  k_constr <- 2
  j_constr <- 1
  p <- 2

  # constraint_fn <- rep(list(function(x){mean(x)}), 2)
  constraint_fn <- rep(list(function(x) pseudohuber_median(x, 0.1)), 2)
  # constraint_grad_fn <- function(x){dpseudohuber_median_dx(x,0.1)
  constraint_grad_fn <- rep(
    list(function(x) {
      dpseudohuber_median_dx(x, 0.1)
    }),
    2
  )

  full_fit <- emuFit_micro_penalized(
    X = X,
    Y = Y,
    B = NULL,
    constraint_fn = constraint_fn,
    tolerance = 1e-7,
    verbose = FALSE
  )

  B <- full_fit$B
  Y_aug <- full_fit$Y_augmented

  X_cup <- X_cup_from_X(X, J)

  j_ref <- 5
  null_fit <- fit_null(
    B = B,
    Y = Y_aug,
    X = X,
    X_cup = X_cup,
    k_constr = k_constr,
    j_constr = j_constr,
    j_ref = j_ref,
    constraint_fn = constraint_fn,
    constraint_grad_fn = constraint_grad_fn,
    constraint_tol = 1e-5,
    B_tol = 1e-4,
    verbose = FALSE,
    trackB = FALSE
  ) ## just track for one j

  null_repar_fit <- fit_null_symmetric(
    B = B,
    Y = Y_aug,
    X = X,
    j_constr = j_constr,
    k_constr = k_constr,
    j_ref = j_ref,
    constraint_fn = constraint_fn,
    constraint_grad_fn = constraint_grad_fn,
    B_tol = 1e-4,
    verbose = FALSE,
    use_optim = TRUE
  )

  null_repar_fit_fs <- fit_null_symmetric(
    B = B,
    Y = Y_aug,
    X = X,
    j_constr = j_constr,
    k_constr = k_constr,
    j_ref = j_ref,
    constraint_fn = constraint_fn,
    constraint_grad_fn = constraint_grad_fn,
    B_tol = 1e-4,
    verbose = FALSE,
    use_optim = FALSE
  )

  #min_mse_lag finds the lagrange multiplier that minimizes the squared norm
  #of the derivative of the lagrangian (ll + lambda*(g - B[k_constr,j_constr])
  #and returns the minimized squared norm -- lower means more accurate fit under null

  null_min_lag_norm <- min_mse_lag(
    X = X,
    Y = Y,
    B = null_fit$B,
    constraint_grad_fn = constraint_grad_fn[[1]],
    k_constr = k_constr,
    j_constr = j_constr,
    j_ref = j_ref
  )
  null_repar_min_lag_norm <- min_mse_lag(
    X = X,
    Y = Y,
    B = null_repar_fit$B,
    constraint_grad_fn = constraint_grad_fn[[1]],
    k_constr = k_constr,
    j_constr = j_constr,
    j_ref = j_ref
  )
  null_repar_min_lag_norm_fs <- min_mse_lag(
    X = X,
    Y = Y,
    B = null_repar_fit_fs$B,
    constraint_grad_fn = constraint_grad_fn[[1]],
    k_constr = k_constr,
    j_constr = j_constr,
    j_ref = j_ref
  )

  #solns are at least close to equal
  # expect_equal(null_repar_fit$B, null_fit$B, tolerance = 1e-2)
  expect_equal(null_repar_fit_fs$B, null_fit$B, tolerance = 1e-2)

  #and to extent that they are not equal, repar fit is more accurate
  # expect_true(null_min_lag_norm > null_repar_min_lag_norm)
  expect_true(null_min_lag_norm > null_repar_min_lag_norm_fs)
})

test_that("compare timing old null fit and symmetric null fit", {
  skip(
    "skip in automated testing because this is a slower bigger simulation study"
  )

  n <- 50
  Js <- c(10, 50, 250)
  ps <- c(2, 4, 8)
  nsim <- 1

  dat <- data.frame(
    bin = rep(0:1, each = 25),
    cont = rnorm(50),
    count = rpois(50, 30) - 30,
    cat = rep(c("A", "B", "C", "D", "E"))
  )
  full_X <- model.matrix(~ bin + cont + count + cat, dat)

  res <- expand.grid(
    seed = 1:nsim,
    J = Js,
    p = ps,
    old_time = NA,
    new_time = NA
  )
  sim_settings <- expand.grid(J = Js, p = ps)

  for (s in 1:nrow(sim_settings)) {
    print(sim_settings[s, ])
    J <- sim_settings$J[s]
    p <- sim_settings$p[s]
    Bs <- get_sim_bs(J)
    B <- rbind(Bs$b0, Bs$b1)
    if (p == 4) {
      B <- rbind(B, rnorm(J), rnorm(J))
      X <- full_X[, 1:4]
    } else if (p == 8) {
      B <- rbind(B, rnorm(J), rnorm(J), rnorm(J), rnorm(J), rnorm(J), rnorm(J))
      X <- full_X
    } else {
      X <- full_X[, 1:2]
    }

    for (sim in 1:nsim) {
      print(sim)
      ind <- which(res$seed == sim & res$J == J & res$p == p)

      Y <- simulate_data(
        n = n,
        J = J,
        distn = "Poisson",
        mean_z = 20,
        B = B,
        X = X
      )
      constraint_fn <- rep(list(function(x) pseudohuber_median(x)), p)
      constraint_grad_fn <- rep(
        list(function(x) {
          dpseudohuber_median_dx(x)
        }),
        p
      )

      full_fit <- try(emuFit_micro_penalized(
        X = X,
        Y = Y,
        B = NULL,
        constraint_fn = constraint_fn,
        #tolerance = 1e-7,
        verbose = FALSE
      ))
      
      if (!inherits(full_fit, "try-error")) {
        B_est <- full_fit$B
        Y_aug <- full_fit$Y_augmented
        
        X_cup <- X_cup_from_X(X, J)
        
        j_ref <- get_j_ref(Y_aug)
        
        print("fitting null")
        start <- proc.time()
        null_fit <- try(fit_null(
          B = B_est,
          Y = Y_aug,
          X = X,
          X_cup = X_cup,
          k_constr = 2,
          j_constr = J / 4,
          j_ref = j_ref,
          constraint_fn = constraint_fn,
          #constraint_tol = 1e-5,
          #B_tol = 1e-4,
          constraint_grad_fn = constraint_grad_fn,
          verbose = FALSE,
          trackB = FALSE
        ))
        end <- proc.time() - start
        if (!inherits(null_fit, "try-error")) {
          res$old_time[ind] <- end[3]
        }
        
        print("fitting new null")
        start <- proc.time()
        null_repar_fit <- try(fit_null_symmetric(
          Y = Y_aug,
          X = X,
          B = B,
          j_constr = J / 4,
          k_constr = 2,
          j_ref = j_ref,
          constraint_fn = constraint_fn[[1]],
          constraint_grad_fn = constraint_grad_fn[[1]],
          #B_tol = 1e-4,
          verbose = TRUE,
          maxit = 1000
        ))
        end <- proc.time() - start
        if (!inherits(null_repar_fit, "try-error")) {
          res$new_time[ind] <- end[3]
        }
      }
    }
  }

  # result here is that timing is pretty similar (double check this)
})

# # profile optim in fit_null_symmetric
# test_that("profile how much time is optim in fit_null_symmetric", {
#   skip("Don't profile new null code when running automatic tests")
# 
#   set.seed(59542234)
#   n <- 100
#   J <- 50
#   X <- cbind(1, rep(c(0, 1), each = n / 2))
#   b0 <- rnorm(J)
#   b1 <- seq(1, 10, length.out = J)
#   b1 <- b1 - mean(b1)
#   b0 <- b0 - mean(b0)
#   Y <- radEmu:::simulate_data(
#     n = n,
#     J = J,
#     X = X,
#     b0 = b0,
#     b1 = b1,
#     distn = "Poisson",
#     mean_z = 8
#   )
# 
#   k_constr <- 2
#   j_constr <- 1
#   p <- 2
# 
#   # constraint_fn <- rep(list(function(x){mean(x)}), 2)
#   constraint_fn <- rep(list(function(x) pseudohuber_median(x, 0.1)), 2)
#   ##### Arguments to fix:
# 
#   # constraint_grad_fn <- function(x){dpseudohuber_median_dx(x,0.1)
#   constraint_grad_fn <- rep(
#     list(function(x) {
#       dpseudohuber_median_dx(x, 0.1)
#     }),
#     2
#   )
# 
#   full_fit <- #suppressMessages(
#     emuFit_micro_penalized(
#       X = X,
#       Y = Y,
#       B = NULL,
#       tolerance = 1e-6,
#       verbose = FALSE
#     )
# 
#   B <- full_fit$B
#   Y_aug <- full_fit$Y_augmented
# 
#   X_cup <- X_cup_from_X(X, J)
# 
#   j_ref <- 5
# 
#   Rprof("../out.prof", interval = 0.01)
#   null_repar_fit <- fit_null_symmetric(
#     Y = Y_aug,
#     X = X,
#     B = B,
#     j_constr = j_constr,
#     k_constr = k_constr,
#     j_ref = j_ref,
#     constraint_fn = constraint_fn[[1]],
#     constraint_grad_fn = constraint_grad_fn[[1]],
#     B_tol = 1e-4,
#     verbose = TRUE,
#     maxit = 1000
#   )
#   Rprof(NULL)
#   summaryRprof("../out.prof")
# })
