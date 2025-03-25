test_that("Penalized estimation reduces to Haldane correction in saturated case", {
  Y <- matrix(1:4,nrow = 2)
  X <- cbind(1,c(0,1))
  penalized_fit <- emuFit_micro_penalized(X,
                                          Y,
                                          B = NULL,
                                          constraint_fn = function(x) x[2],
                                          maxit = 500,
                                          tolerance = 1e-5,
                                          verbose = FALSE)
  
  ml_fit <- emuFit_micro(X,Y + 0.5,
                         constraint_fn = function(x) x[2],
                         tolerance = 1e-5,
                         verbose = FALSE)
  
  expect_true(max(abs(penalized_fit$B- ml_fit)) <1e-4)
  expect_equal(penalized_fit$Y_augmented - Y, matrix(0.5, nrow = 2, ncol = 2))
  
  
  
  
  
})

test_that("Penalized estimation reduces to Haldane correction in saturated case with zeroes", {
  Y <- matrix(0:3,nrow = 2)
  X <- cbind(1,c(0,1))
  penalized_fit <- emuFit_micro_penalized(X,
                                          Y,
                                          B = NULL,
                                          constraint_fn = function(x) x[2],
                                          maxit = 500,
                                          tolerance = 1e-5,
                                          verbose = FALSE)
  
  ml_fit <- emuFit_micro(X,Y + 0.5,
                         constraint_fn = function(x) x[2],
                         tolerance = 1e-5,
                         verbose = FALSE)
  
  expect_true(max(abs(penalized_fit$B - ml_fit))<1e-4)
  expect_equal(penalized_fit$Y_augmented - Y, matrix(0.5, nrow = 2, ncol = 2))
})



test_that("PL fit to simple example returns reasonable values", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  J <- 10
  n <- 40
  b1 <- seq(1,5,length.out = J) -
    mean(seq(1,5,length.out = J))
  Y <- radEmu:::simulate_data(n = n,
                              J = J,
                              X = X,
                              b0 = rnorm(J),
                              b1 = b1,
                              distn = "Poisson",
                              mean_z = 50)
  
  pl_fit <- emuFit_micro_penalized(X,
                                   Y,
                                   B = matrix(rnorm(20),nrow = 2),
                                   constraint_fn = function(x) mean(x),
                                   maxit = 200,
                                   tolerance = 1e-10,
                                   verbose= FALSE)
  
  # plot(b1-mean(b1),ml_fit[2,])
  # abline(a = 0,b = 1,lty =2,col = "red")
  
  expect_true(max(abs(pl_fit$B[2,] - (b1 - mean(b1))))<1e-2)
})


test_that("PL fit to simple example returns numerically identical results
regardless of whether we use computationally efficient augmentation or older
less efficient implementation (and that both substantially differ from MLE", {
  set.seed(4323)
  J <- 10
  n <- 10
  X <- cbind(1,rep(c(0,1),each = n/2))
  Y <- radEmu:::simulate_data(n = n,
                              J = J,
                              X = X,
                              b0 = rnorm(J),
                              b1 = seq(1,5,length.out = J),
                              distn = "ZINB",
                              zinb_size = 2,
                              zinb_zero_prop = 0.7,
                              mean_z = 10)
  
  ml_fit <-  emuFit_micro(X = as.matrix(X),
                          Y = as.matrix(Y),
                          B = NULL,
                          constraint_fn = function(x) mean(x),
                          maxit = 1000,
                          tolerance = 1e-3,
                          verbose= FALSE)
  
  pl_fit_new <- emuFit_micro_penalized(X,
                                       Y,
                                       B = NULL,
                                       constraint_fn = function(x) mean(x),
                                       maxit = 1000,
                                       tolerance = 1e-3,
                                       verbose= FALSE)
  #using old implementation should trigger a message explaining there's no
  #reason to do this except for testing
  suppressMessages(expect_message(
    pl_fit_old <- emuFit_micro_penalized(X,
                                         Y,
                                         B = NULL,
                                         constraint_fn = function(x) mean(x),
                                         maxit = 10000,
                                         tolerance = 1e-3,
                                         verbose= TRUE,
                                         use_legacy_augmentation = TRUE)))
  
  expect_true(max(abs(ml_fit - pl_fit_new$B))>1)
  expect_true(max(abs(ml_fit - pl_fit_old$B))>1)
  expect_equal(pl_fit_new,pl_fit_old,tolerance = 1e-6)
})



test_that("PL fit to simple continuous covariate example returns numerically identical results
regardless of whether we use computationally efficient augmentation or older
less efficient implementation (and that both substantially differ from MLE", {
  set.seed(4323)
  X <- cbind(1,rnorm(10))
  J <- 10
  n <- 10
  Y <- radEmu:::simulate_data(n = n,
                              J = J,
                              X = X,
                              b0 = rnorm(J),
                              b1 = seq(1,5,length.out = J),
                              distn = "ZINB",
                              zinb_size = 2,
                              zinb_zero_prop = 0.7,
                              mean_z = 10)
  
  ml_fit <-  emuFit_micro(X,
                          Y,
                          B = matrix(0,nrow = 2, ncol = J),
                          constraint_fn = function(x) mean(x),
                          maxit = 10000,
                          tolerance = 0.01,
                          verbose= FALSE)
  
  pl_fit_new <- emuFit_micro_penalized(X,
                                       Y,
                                       B = NULL,
                                       constraint_fn = function(x) mean(x),
                                       maxit = 10000,
                                       tolerance = 0.001,
                                       verbose= FALSE)
  #using old implementation should trigger a message explaining there's no
  #reason to do this except for testing
  suppressMessages(expect_message(
    pl_fit_old <- emuFit_micro_penalized(X,
                                         Y,
                                         B = NULL,
                                         constraint_fn = function(x) mean(x),
                                         maxit = 10000,
                                         tolerance = 0.001,
                                         verbose= TRUE,
                                         use_legacy_augmentation = TRUE)))
  
  expect_true(max(abs(ml_fit - pl_fit_new$B))>0.5)
  expect_true(max(abs(ml_fit - pl_fit_old$B))>0.5)
  expect_equal(pl_fit_new,pl_fit_old)
})

test_that("penalized fit uses B if given, and therefore fit is quicker", {
  set.seed(4323)
  X <- cbind(1,rnorm(10))
  J <- 10
  n <- 10
  Y <- radEmu:::simulate_data(n = n,
                              J = J,
                              X = X,
                              b0 = rnorm(J),
                              b1 = seq(1,5,length.out = J),
                              distn = "ZINB",
                              zinb_size = 2,
                              zinb_zero_prop = 0.7,
                              mean_z = 10)
  
  start <- proc.time()
  pl_fit_one <- emuFit_micro_penalized(X,
                                       Y,
                                       B = NULL,
                                       constraint_fn = function(x) mean(x),
                                       maxit = 10000,
                                       tolerance = 0.01,
                                       verbose= FALSE)
  end <- proc.time() - start 
  start_refit <- proc.time() 
  pl_fit_two <- emuFit_micro_penalized(X,
                                       Y,
                                       B = pl_fit_one$B,
                                       constraint_fn = function(x) mean(x),
                                       maxit = 10000,
                                       tolerance = 0.01,
                                       verbose= FALSE)
  end_refit <- proc.time() - start_refit 
  expect_true(end_refit[3] < end[3])
  
})




