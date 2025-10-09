test_that("in saturated case, augmentations reduce to haldane correction", {
    Y <- matrix(1:4,nrow = 2)
    X <- cbind(1,c(0,1))
    J <- ncol(Y)
    n <- nrow(Y)


    ml_fit <- emuFit_micro(X,Y ,
                           constraint_fn = rep(list(function(x) x[2]), 2), 
                           tolerance = 1e-5,
                           verbose = FALSE)

    X_cup <- X_cup_from_X_fast(X,J)
    G <- get_G_for_augmentations_fast(X = X,J = J,n = n,
                                 X_cup = X_cup)
    expect_true(max(abs(get_augmentations(X = X,
                                          G = G,
                                          Y = Y,
                                          B = ml_fit) - 0.5)) < 1e-10)

    #correction should not depend on B in this case -- check this too
    expect_true(max(abs(get_augmentations(X = X,
                                          G = G,
                                          Y = Y,
                                          B = 0*ml_fit) - 0.5)) < 1e-10)

  })

test_that("get_augmentations_par gives same results as get_augmentations", {
  
  set.seed(59542234)
  n <- 50
  J <- 10
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
  
  X_cup <- X_cup_from_X_fast(X, ncol(Y))
  G <- get_G_for_augmentations_fast(X, ncol(Y), nrow(Y), X_cup)
  j_ref <- get_j_ref(Y)
  Y_augmented <- as.matrix(Y) + 1e-3*mean(as.numeric(unlist(Y)))
  fitted_model <- emuFit_micro(X,
                               Y_augmented,
                               j_ref = j_ref, maxit = 5)
  aug1 <- get_augmentations(X, G, Y, fitted_model)
  aug2 <- get_augmentations_par(X, G, Y, fitted_model, par = FALSE)
  aug3 <- get_augmentations_par(X, G, Y, fitted_model, par = TRUE)

  expect_true(all.equal(aug1, aug2))
  expect_true(all.equal(aug2, aug3))
  
})

test_that("get_augmentations_par is faster for large n", {
  
  skip("test takes too long to be automated")
  
  set.seed(59542234)
  n <- 1000
  J <- 100
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
  
  X_cup <- X_cup_from_X_fast(X, ncol(Y))
  G <- get_G_for_augmentations_fast(X, ncol(Y), nrow(Y), X_cup)
  j_ref <- get_j_ref(Y)
  Y_augmented <- as.matrix(Y) + 1e-3*mean(as.numeric(unlist(Y)))
  fitted_model <- emuFit_micro(X,
                               Y_augmented,
                               j_ref = j_ref, maxit = 5)
  start1 <- proc.time()
  aug1 <- get_augmentations(X, G, Y, fitted_model)
  end1 <- proc.time() - start1
  
  start2 <- proc.time()
  aug2 <- get_augmentations_par(X, G, Y, fitted_model, par = FALSE)
  end2 <- proc.time() - start2
  
  start3 <- proc.time()
  aug3 <- get_augmentations_par(X, G, Y, fitted_model, par = TRUE)
  end3 <- proc.time() - start3
  
  if (.Platform$OS.type == "unix" & parallel::detectCores() > 4) {
    expect_true(end3[3] < end1[3])
    expect_true(end3[3] < end2[3])
  }
  
})
