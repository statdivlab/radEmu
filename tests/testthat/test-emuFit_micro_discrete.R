
test_that("debug first implementation", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  J <- 10
  n <- 40
  b1 <- 1:J - mean(1:J)
  Y <- radEmu:::simulate_data(n = n,
                              J = J,
                              X = X,
                              b0 = rnorm(J),
                              b1 = b1,
                              distn = "Poisson",
                              mean_z = 10)
  
  ml_fit <- emuFit_micro(X = X,
                         Y = Y,
                         constraint_fn = rep(list(function(x) mean(x)), 2), 
                         maxit = 200,
                         tolerance = 1e-3,
                         verbose = FALSE)
  ml_fit_fisher <- emuFit_micro(X = X,
                                Y = Y,
                                constraint_fn = rep(list(function(x) mean(x)), 2), 
                                maxit = 200,
                                tolerance = 1e-3,
                                verbose = FALSE,
                                use_discrete = FALSE)
  
  
  amys_disaster <- emuFit_micro_discrete(X = X,
                                         Y = Y,
                                         j_ref = 1)
  
  
  testthat::expect_lte(  max(abs(ml_fit -  
                                    rbind(amys_disaster[1,] - mean(amys_disaster[1,]), 
                                          amys_disaster[2,] - mean(amys_disaster[2,])))),
                          5e-3)
  
  testthat::expect_lte(max(abs(ml_fit - ml_fit_fisher)), 5e-3)
  
})
