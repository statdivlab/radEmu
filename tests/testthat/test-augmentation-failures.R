test_that("confirm Matrix Csparse_transpose issue is not happening", {
  
  set.seed(1)
  X <- structure(c(rep(1, 18), rep(0, 6)), dim = c(12L, 2L))
  
  Y <- radEmu:::simulate_data(n = 12L, J = 4L,
                              b0 = runif(4L, min = 0, max = 4),
                              b1 = runif(4L, min = 0, max = 4),
                              X = X,
                              distn = "Poisson",
                              mean_z = 5)
  
  covariates <- structure(list(group = c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)),
                          class = "data.frame", row.names = c(NA, -12L))
  
  # devtools::load_all()  # info <- methods::as(info, "symmetricMatrix")
  fitted_model <- emuFit(Y = Y,
                         X = X,
                         formula = ~group,
                         data = covariates,
                         verbose = FALSE,
                         B_null_tol = 1e-2,
                         tolerance = 0.01,
                         tau = 2,
                         run_score_test = TRUE,
                         return_wald_p = TRUE)
  
  expect_true("emuFit" %in% class(fitted_model))
  
  
  ### check data frame inputs ok
  fitted_model_df <- emuFit(Y = as.data.frame(Y)[5:8, ],
                            X = as.data.frame(X)[5:8, ],
                            formula = ~group,
                            data = covariates,
                            verbose = FALSE,
                            B_null_tol = 1e-2,
                            tolerance = 0.01,
                            tau = 2,
                            run_score_test = TRUE,
                            return_wald_p = TRUE)
  
  expect_true("emuFit" %in% class(fitted_model_df))
  
})
