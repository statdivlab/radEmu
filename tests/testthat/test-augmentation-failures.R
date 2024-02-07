test_that("confirm Matrix Csparse_transpose issue is not happening", {
  
  Y <- structure(c(1087, 3541, 0, 2432, 0, 18538, 1158, 2282, 625, 0, 
                   3759, 0, 0, 4658, 0, 3719, 0, 0, 7531, 48316, 0, 0, 20273, 0, 
                   1227, 0, 1471, 0, 479, 10602, 3115, 3286, 0, 1969, 0, 3045, 0, 
                   8018, 0, 1622, 1307, 34117, 9338, 0, 0, 7909, 0, 0), dim = c(12L, 4L))
  X <- structure(c(rep(1, 18), rep(0, 6)), dim = c(12L, 2L))
  covariates <- structure(list(group = c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)), class = "data.frame", row.names = c(NA, -12L))
  
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
