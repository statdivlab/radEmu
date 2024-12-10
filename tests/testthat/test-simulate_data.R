test_that("simulating data works", {
  J <- 10; n <- 60
  # cluster membership 
  cluster <- rep(1:4, each = n/4)
  # intercepts for each category
  b0 <- rnorm(J)
  # coefficients for X1 for each category 
  b1 <- seq(1, 5, length.out = J)
  
  dat1 <- radEmu:::simulate_data(n = n, J = J, b0 = b0, b1 = b1, distn = "Poisson", 
                                 mean_z = 8)
  dat2 <- radEmu:::simulate_data(n = n, J = J, b0 = b0, b1 = b1, distn = "Poisson", 
                                 mean_z = 8, cluster = cluster)
  dat3 <- radEmu:::simulate_data(n = n, J = J, b0 = b0, b1 = b1, distn = "ZINB", 
                                 mean_z = 8, zinb_size = 10, zinb_zero_prop = 0.5)
  dat4 <- radEmu:::simulate_data(n = n, J = J, b0 = b0, b1 = b1, distn = "ZINB", cluster = cluster,
                                 mean_z = 8, zinb_size = 10, zinb_zero_prop = 0.5)
  
  # make sure data are the correct dimensions
  expect_true(nrow(dat1) == n & nrow(dat2) == n & nrow(dat3) == n & nrow(dat4) == n &
                ncol(dat1) == J & ncol(dat2) == J & ncol(dat3) == J & ncol(dat4) == J)
  # check that ZINB data is more sparse than Poisson data
  expect_true(mean(dat3 == 0) > mean(dat1 == 0))
  expect_true(mean(dat4 == 0) > mean(dat2 == 0))
  
})

test_that("simulating b's works", {
  b_res <- get_sim_bs(10)
  expect_true(length(b_res[[1]]) == 10) 
  expect_true(length(b_res[[2]]) == 10)
})
  