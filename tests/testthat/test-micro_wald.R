
# Here, we manually insert data used in test before the simulate_data() function
# was changed. This is a temporary solution until we want to do away with the
# specific p-value confirmations
Y_old1 <-
  matrix(
    c( 1748,      0,     0,  1287,     0,     0,      0,      0,   1875,     0,
       8286, 124134,     0,  1716,     0,     0,   4529,   3564,  21975,     0,
       4096,  43569,     0,     0,     0,     0,   2885,   2250,   1685,   217,
       4289, 122134,   416,     0,  5936,     0,      0,      0,      0,     0,
       1122,  15785,   920,  1354,     0,     0,      0,      0,   1654,     0,
       30007,  99540,     0,  8783,     0,     0,  13382,  20486,  28670,     0,
       5087,      0,     0,  3040,     0,  3954,  11802,      0,   9331,     0,
       3841,  41104,  1931,  6271, 33557,  1338,      0,   3442,      0,     0,
       3059,      0,  1279,  2274, 11459,   886,   2144,   5133,   4765,     0,
       3105,      0,     0,     0,     0,   426,   2622,   5103,      0,     0, 
         80,      0,     0, 26431,     0,     0,  92214, 299927,      0,     0, 
         32,      0,     0,  5186,  4065,     0,  18326,  34273,      0, 76950,
          0,      0,     2,  4147,     0,     0,   6183,  17407,  76158,     0,
         20,    572,    21,     0,  5391,   496,  11737,      0,      0, 19911,
         13,      0,    49,     0,  6721,     0,      0,  17896,  85495, 33042,
          0,   1497,    41,  6450,  8997,   709,      0, 149402, 247396, 43690,
          0,      0,     0,     0,  9225,   508,  12808,  42592,  51659, 68014,
         30,      0,    89,     0, 13951,   840,   7604,      0,  93587, 43885,
         41,    314,     0,  1483,  4061,   680,   4684,      0,      0,  6086,
         54,      0,    85,     0,  3871,     0,   9348,  25395,  84277,     0),
    nrow = 20, ncol = 10, byrow = TRUE)

Y_old2 <-
  matrix(
    c(6051354, 628305, 417132,     0,    0,  20886,     0,    237,      0,     54,
       190080,      0,      0,  327,  234,   5265,   106,    177,      0,     45,
            0,      0,  37588,  578,  933,  16147,     0,    691,    906,   1277,
            0,      0,      0,     0,    0,   1788,  1122,      0,   7657,  34478,
            0,      0,  16470,  9412,    0,  23456,     0,  17289,      0,      0,
            0,   1459,      0,  1479, 3321,  14736,  1319,      0,      0,      0,
           13,     13,      0,    324, 1832,      0,     0,      0,      0,      0,
            0,  22292,  21209,     0,    0,   2387,   139,    109,      6,      0,
            0,  17457,      0,  715, 1298,      0,   245,      0,      0,    127,
            0,   2562,      0,     0,    0,   1752,   522,      0,      0,      0,
          153,    113,    226,     0,  566,   2176,   861,      0,      0,      0,
            0,   6512, 110014, 10616, 6482,  86103,  2774,  18807,   7569,   8734,
          457,      0,   1240,   351,  401,      0,  3035,  10084,   4080,  14265,
      2631304, 255618,      0, 2257,    0,      0,  1351,    798,      0,      0,
            0,      0,      0,     0,    0,      0,   145,      0,      5,      0,
       110612,   8278,      0,  1309,    0,   5002,     0,   1018,      0,      0,
            0,      0,      0,     0, 5779, 230361,     0, 431527, 450787,      0,
            0,      0,  22245,   139,    0,      0,   151,      0,      0,      0,
          180,    169,   2659,   1115,    0,   6693, 12595, 125758, 107518, 425831,
            0,      0,  20670,  1216, 2086,  13617,  3777,      0,   1562,   4535),
    nrow = 20, ncol = 10, byrow = TRUE)


test_that("wald test gives semi-reasonable output with categorical covariate", {
  set.seed(343234)
  n <- 20
  X <- cbind(1,rep(c(0,1),each = n/2))
  J <- 10
  b0 <- rnorm(10)
  b1 <- 1:10
  b1 <- b1 - mean(b1)
  b1[5] <- pseudohuber_center(b1[-5],0.1)
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
                              mean_z = 8)
  
  k_constr <- 2
  j_constr <- 5
  p <- 2
  
  constraint_fn <- function(x){ pseudohuber_center(x,0.1)}
  
  ##### Arguments to fix:
  
  constraint_grad_fn <- function(x){dpseudohuber_center_dx(x,0.1)}
  
  constraint_hess_fn <- function(x,ind_1,ind_2){hess_pseudohuber_center(x,0.1,ind_1,ind_2)}
  b[2,4] <- constraint_fn(b[2,-4])
  
  X_cup <- X_cup_from_X(X,J)
  
  full_fit <- emuFit_micro_penalized(X = X,
                                     Y = Y_old1,
                                     B = NULL,
                                     constraint_fn = constraint_fn,
                                     tolerance = 1e-3,
                                     verbose = FALSE)
  
  wald_result <-
    micro_wald(Y = full_fit$Y,
               X,
               X_cup = X_cup,
               B = full_fit$B,
               test_kj = data.frame(k = 2, j = 4),
               constraint_fn = constraint_fn,
               constraint_grad_fn = constraint_grad_fn,
               nominal_coverage = 0.95)
  
  expect_true(is.data.frame(wald_result$coefficients))
  expect_true(wald_result$coefficients$pval>0.1)
  expect_true(is.list(wald_result))
  expect_true(ncol(wald_result$I) ==20)
  expect_equal(wald_result$coefficients$pval, 0.61, tolerance = 0.02)
  
  
  
  full_fit_new <- emuFit_micro_penalized(X = X,
                                         Y = Y,
                                         B = NULL,
                                         constraint_fn = constraint_fn,
                                         tolerance = 1e-3,
                                         verbose = FALSE)
  
  wald_result_new <-
    micro_wald(Y = full_fit_new$Y,
               X,
               X_cup = X_cup,
               B = full_fit_new$B,
               test_kj = data.frame(k = 2, j = 4),
               constraint_fn = constraint_fn,
               constraint_grad_fn = constraint_grad_fn,
               nominal_coverage = 0.95)
  
  expect_true(is.data.frame(wald_result_new$coefficients))
  expect_true(wald_result_new$coefficients$pval>0.1)
  expect_true(is.list(wald_result_new))
  expect_true(ncol(wald_result_new$I) ==20)
  
})



test_that("wald test gives semi-reasonable output with continuous covariate", {
  set.seed(9944234)
  n <- 20
  X <- cbind(1,rnorm(n))
  J <- 10
  b0 <- rnorm(10)
  b1 <- 1:10
  b1 <- b1 - mean(b1)
  b1[5] <- pseudohuber_center(b1[-5],0.1)
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
                              mean_z = 8)
  
  k_constr <- 2
  j_constr <- 5
  p <- 2
  
  constraint_fn <- function(x){ pseudohuber_center(x,0.1)}
  
  ##### Arguments to fix:
  
  constraint_grad_fn <- function(x){dpseudohuber_center_dx(x,0.1)}
  
  constraint_hess_fn <- function(x,ind_1,ind_2){hess_pseudohuber_center(x,0.1,ind_1,ind_2)}
  b[2,4] <- constraint_fn(b[2,-4])
  
  X_cup <- X_cup_from_X(X,J)
  
  full_fit <- emuFit_micro_penalized(X = X,
                                     Y = Y_old2,
                                     B = NULL,
                                     constraint_fn = constraint_fn,
                                     tolerance = 1e-3,
                                     verbose = FALSE)
  
  wald_result <- micro_wald(Y = full_fit$Y,
                            X,
                            X_cup = X_cup,
                            B = full_fit$B,
                            test_kj = data.frame(k = 2, j = 4),
                            constraint_fn = constraint_fn,
                            constraint_grad_fn = constraint_grad_fn,
                            nominal_coverage = 0.95)
  
  wald_result_for_an_alternative <- micro_wald(Y = full_fit$Y,
                                               X,
                                               X_cup = X_cup,
                                               B = full_fit$B,
                                               test_kj = data.frame(k = 2, j = 10),
                                               constraint_fn = constraint_fn,
                                               constraint_grad_fn = constraint_grad_fn,
                                               nominal_coverage = 0.95)
  
  
  expect_true(is.data.frame(wald_result$coefficients))
  expect_true(is.list(wald_result))
  expect_true(ncol(wald_result$I) ==20)
  expect_equal(wald_result$coefficients$pval, 0.11, tolerance = 0.03)
  expect_true(wald_result_for_an_alternative$coefficients$pval < 0.01)
  
  
  
  full_fit_new <- emuFit_micro_penalized(X = X,
                                     Y = Y,
                                     B = NULL,
                                     constraint_fn = constraint_fn,
                                     tolerance = 1e-3,
                                     verbose = FALSE)
  
  wald_result_new <- micro_wald(Y = full_fit_new$Y,
                                X,
                                X_cup = X_cup,
                                B = full_fit_new$B,
                                test_kj = data.frame(k = 2, j = 4),
                                constraint_fn = constraint_fn,
                                constraint_grad_fn = constraint_grad_fn,
                                nominal_coverage = 0.95)
      
  wald_result_for_an_alternative_new <- micro_wald(Y = full_fit_new$Y,
                                                   X,
                                                   X_cup = X_cup,
                                                   B = full_fit_new$B,
                                                   test_kj = data.frame(k = 2, j = 10),
                                                   constraint_fn = constraint_fn,
                                                   constraint_grad_fn = constraint_grad_fn,
                                                   nominal_coverage = 0.95)
  
  
  expect_true(is.data.frame(wald_result_new$coefficients))
  expect_true(is.list(wald_result_new))
  expect_true(ncol(wald_result_new$I) ==20)
  expect_true(wald_result_for_an_alternative_new$coefficients$pval < 0.01)
  
})
