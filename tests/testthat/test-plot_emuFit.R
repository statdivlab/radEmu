set.seed(11)
J <- 6
n <- 12
X <- cbind(1,rnorm(n))
Y <- radEmu:::simulate_data(n = n,
                            J = J,
                            X = X,
                            b0 = rnorm(10),
                            b1 = 1:10 - mean(1:10),
                            distn = "ZINB",
                            zinb_size = 3,
                            zinb_zero_prop = 0.6,
                            mean_z = 8)

fitted_model <- emuFit(Y = Y,
                       X = X,
                       verbose = FALSE,
                       B_null_tol = 1e-2,
                       tolerance = 0.01,
                       tau = 2,
                       return_wald_p = FALSE,
                       compute_cis = TRUE,
                       run_score_tests = FALSE, 
                       use_fullmodel_info = FALSE,
                       use_fullmodel_cov = FALSE,
                       return_both_score_pvals = FALSE)

test_that("plot() returns data frame and plot", {
  plot_out <- plot(x = fitted_model)
  expect_true(is.data.frame(plot_out$data))
  expect_true(all(sapply(plot_out$plots, ggplot2::is.ggplot)))
})


test_that("plot() returns error when plot_key does not match coefficient table", {
  expect_error({
    plot(x = fitted_model,
         plot_key = list(c("First Covariate" = "covariate1")))
  })
})

test_that("plot() returns error when coefficient is included multiple times in plot_key", {
  expect_error({
    plot(x = fitted_model,
         plot_key = list(c("First Covariate" = "covariate_1"),
                         c("Second Covariate" = "covariate_1")))
  })
})



