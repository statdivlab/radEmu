set.seed(11)
J <- 6
p <- 2
n <- 12
X <- cbind(1,rnorm(n))
z <- rnorm(n) +5
b0 <- rnorm(J)
b1 <- seq(1,5,length.out = J)
b1 <- b1 - mean(b1)
b <- rbind(b0,b1)
Y <- matrix(NA,ncol = J, nrow = n)

for(i in 1:n){
  for(j in 1:J){
    temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
    Y[i,j] <- rnbinom(1, mu= temp_mean,size = 2)*rbinom(1,1,0.8)
  }
}

fitted_model <- emuFit(Y = Y,
                       X = X,
                       formula = ~group,
                       data = covariates,
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

test_that("plot.radEmu returns data frame and plot", {
  plot_out <- plot.radEmu(x = fitted_model)
  expect_true(is.data.frame(plot_out$data))
  expect_true(all(sapply(plot_out$plots, ggplot2::is.ggplot)))
})


test_that("plot.radEmu returns error when plot_key does not match coefficient table", {
  expect_error({
    plot.radEmu(x = fitted_model,
                plot_key = list(c("First Covariate" = "covariate1")))
  })
})

test_that("plot.radEmu returns error when coefficient is included multiple times in plot_key", {
  expect_error({
    plot.radEmu(x = fitted_model,
                plot_key = list(c("First Covariate" = "covariate_1"),
                                c("Second Covariate" = "covariate_1")))
  })
})



