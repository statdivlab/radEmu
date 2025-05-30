
test_that("ML fit to simple example give reasonable output", {
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

  # plot(b1-mean(b1),ml_fit[2,])
  # abline(a = 0,b = 1,lty =2,col = "red")

  expect_true(max(abs(ml_fit[2,] - (b1 - mean(b1)))) < 1) # increased to 1
})

test_that("With or without 'working_constraint' we get same results", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  J <- 10
  n <- 40
  Y <- radEmu:::simulate_data(n = n,
                              J = J,
                              X = X,
                              b0 = rnorm(J),
                              b1 = seq(1,10,length.out = J),
                              distn = "Poisson",
                              mean_z = 5)
  
  ml_fit <- emuFit_micro(X,
                         Y,
                         constraint_fn = rep(list(function(x) mean(x)), 2),
                         maxit = 200,
                         tolerance = 1e-6,
                         verbose= FALSE)
  
  ml_fit_direct <- emuFit_micro(X,
                                Y,
                                constraint_fn = rep(list(function(x) mean(x)), 2),
                                maxit = 200,
                                use_working_constraint = FALSE,
                                tolerance = 1e-6,
                                verbose = FALSE)


  expect_equal(ml_fit,ml_fit_direct,tolerance = 1e-4)
})



test_that("PL fit with categorical predictor matches analytical form of MPLE in this case,
          and does NOT match MLE when group sizes are equal", {
            set.seed(90333)
            X <- cbind(1,rep(c(0,1),each = 20))
            z <- rnorm(40) 
            J <- 10
            p <- 2
            n <- 40
            b0 <- rnorm(J)
            b1 <- seq(1,10,length.out = J)/5
            b <- rbind(b0,b1)
            Y <- matrix(NA,ncol = J, nrow = 40)
            
            for(i in 1:40){
              for(j in 1:J){
                temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
                Y[i,j] <- rpois(1, lambda = temp_mean)
              }
            }
            pl_fit <- emuFit_micro_penalized(X,
                                             Y,
                                             B = matrix(rnorm(20),nrow = 2),
                                             constraint_fn = rep(list(function(x) mean(x)), 2),
                                             maxit = 200,
                                             tolerance = 1e-8,
                                             verbose= FALSE)
            
            ml_fit <- emuFit_micro(X,
                                   Y,
                                   B = matrix(rnorm(20),nrow = 2),
                                   constraint_fn = rep(list(function(x) mean(x)), 2),
                                   maxit = 200,
                                   tolerance = 1e-8,
                                   verbose= FALSE)
            
            cs_grp1 <- colSums(Y[1:20,])
            cs_grp2 <- colSums(Y[21:40,])
            
            bhat1 <- log(cs_grp1 + 0.5) - log(cs_grp1[1] + 0.5)
            bhat2 <- log(cs_grp2 + 0.5) - log(cs_grp2[1] +0.5)
            bhat2 <- bhat2 - bhat1
            bhat1 <- bhat1 - mean(bhat1)
            bhat2 <- bhat2 - mean(bhat2)
            analytical_B <- rbind(bhat1,bhat2)
            
            expect_true(max(abs(pl_fit$B - analytical_B))< 1e-7)
            expect_true(max(abs(ml_fit - analytical_B))>0.01)
          })


test_that("PL fit with categorical predictor matches analytical form of MPLE in this case,
          and does NOT match MLE when group sizes are unequal", {
            set.seed(90333)
            X <- cbind(1,rep(c(0,0,0,1),each = 10))
            z <- rnorm(40) 
            J <- 10
            p <- 2
            n <- 40
            b0 <- rnorm(J)
            b1 <- seq(1,10,length.out = J)/5
            b <- rbind(b0,b1)
            Y <- matrix(NA,ncol = J, nrow = 40)
            
            for(i in 1:40){
              for(j in 1:J){
                temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
                Y[i,j] <- rpois(1, lambda = temp_mean)
              }
            }
            pl_fit <- emuFit_micro_penalized(X,
                                             Y,
                                             B = matrix(rnorm(20),nrow = 2),
                                             constraint_fn = rep(list(function(x) mean(x)), 2),
                                             maxit = 200,
                                             tolerance = 1e-8,
                                             verbose= FALSE)
            
            ml_fit <- emuFit_micro(X,
                                   Y,
                                   B = matrix(rnorm(20),nrow = 2),
                                   constraint_fn = rep(list(function(x) mean(x)), 2),
                                   maxit = 200,
                                   tolerance = 1e-8,
                                   verbose= FALSE)
            
            cs_grp1 <- colSums(Y[1:30,])
            cs_grp2 <- colSums(Y[31:40,])
            
            bhat1 <- log(cs_grp1 + 0.5) - log(cs_grp1[1] + 0.5)
            bhat2 <- log(cs_grp2 + 0.5) - log(cs_grp2[1] +0.5)
            bhat2 <- bhat2 - bhat1
            bhat1 <- bhat1 - mean(bhat1)
            bhat2 <- bhat2 - mean(bhat2)
            analytical_B <- rbind(bhat1,bhat2)
            
            expect_true(max(abs(pl_fit$B - analytical_B))< 1e-7)
            expect_true(max(abs(ml_fit - analytical_B))>0.01)
          })

test_that("We get same results with and without warm start", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  J <- 10
  n <- 40
  Y <- radEmu:::simulate_data(n = n,
                              J = J,
                              X = X,
                              b0 = rnorm(J),
                              b1 = 1:J,
                              distn = "Poisson",
                              mean_z = 2)

  # may need to have large number of iterations and small tolerance
  ml_fit <- emuFit_micro(X,
                         Y,
                         constraint_fn = rep(list(function(x) mean(x)), 2),
                         maxit = 1e3,
                         tolerance = 1e-14,
                         verbose = FALSE)
  
  ml_fit_direct <- emuFit_micro(X,
                                Y,
                                constraint_fn = rep(list(function(x) mean(x)), 2),
                                maxit = 1e3,
                                warm_start = FALSE,
                                tolerance = 1e-14,
                                verbose = FALSE)

  expect_equal(ml_fit, ml_fit_direct, tolerance = 1e-6)
})



test_that("We get a fit if we don't specify constraint", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  J <- 10
  n <- 40
  Y <- radEmu:::simulate_data(n = n,
                              J = J,
                              X = X,
                              b0 = rnorm(J),
                              b1 = 1:J,
                              distn = "Poisson",
                              mean_z = 2)
  
  ml_fit <- emuFit_micro(X,
                         Y,
                         maxit = 200,
                         tolerance = 1e-6,
                         verbose = FALSE)


  expect_true(is.matrix(ml_fit))
})

test_that("ML fit to simple example give reasonable output with J >> n", {
  set.seed(4323)
  n <- 10
  X <- cbind(1,rep(c(0,1),each = n/2))
  J <- 1000
  b0 <- rnorm(J)
  b1 <- seq(-5,5,length.out = J)
  b <- rbind(b0,b1)
  
  Y <- radEmu:::simulate_data(n = n,
                              J = J,
                              X = X,
                              b0 = b0,
                              b1 = b1,
                              distn = "Poisson",
                              mean_z = 10)

  ml_fit <- emuFit_micro(X,
                         Y,
                         constraint_fn = rep(list(function(x) mean(x)), 2),
                         maxit = 500,
                         tolerance = 1e-3,
                         verbose = FALSE)
  
  expect_true(max(abs(ml_fit[2,] - b1))<.5)
})

test_that("unpenalized fit uses B if given, and therefore fit is quicker", {
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
  pl_fit_one <- emuFit_micro(X,
                             Y,
                             B = NULL,
                             constraint_fn = rep(list(function(x) mean(x)), 2),
                             maxit = 10000,
                             tolerance = 0.01,
                             verbose= FALSE)
  end <- proc.time() - start 
  start_refit <- proc.time() 
  pl_fit_two <- emuFit_micro(X,
                             Y,
                             B = pl_fit_one,
                             constraint_fn = rep(list(function(x) mean(x)), 2),
                             maxit = 10000,
                             tolerance = 0.01,
                             verbose= FALSE)
  end_refit <- proc.time() - start_refit 
  expect_true(end_refit[3] < end[3])
  
})

test_that("unpenalized fit converges quicker if optimize_rows is set to TRUE", {
  set.seed(4323)
  J <- 100
  n <- 100
  X <- cbind(1,rnorm(n))
  Y <- radEmu:::simulate_data(n = n,
                              J = J,
                              X = X,
                              b0 = rnorm(J),
                              b1 = seq(1,5,length.out = J),
                              distn = "Poisson",
                              # zinb_size = 2,
                              # zinb_zero_prop = 0.7,
                              mean_z = 10)
  
  start <- proc.time()
  pl_fit_one <- emuFit_micro(X,
                             Y,
                             B = NULL,
                             # constraint_fn = function(x) mean(x),
                             maxit = 10000,
                             tolerance = 1e-6,
                             verbose= FALSE,
                             optimize_rows = FALSE)
  end <- proc.time() - start 
  start_refit <- proc.time() 
  pl_fit_two <- emuFit_micro(X,
                             Y,
                             B = NULL,
                             # constraint_fn = function(x) mean(x),
                             maxit = 10000,
                             tolerance =1e-6,
                             verbose= FALSE,
                             optimize_rows = TRUE)
  end_refit <- proc.time() - start_refit 
  # confirm that new approach is faster
  expect_true(end_refit[3] < end[3])
  # confirm that two estimates are very similar
  expect_true(max(pl_fit_one - pl_fit_two) < 1e-5)
  
})

test_that("unpenalized fit converges quicker if optimize_rows is set to TRUE with large p", {
  
  skip(message = "Skipping because test is very slow with J = 100 and p = 8")
  
  set.seed(4323)
  J <- 100
  n <- 100
  X <- cbind(1,rnorm(n),rnorm(n), rep(c(0, 1, 0, 0), each = 25), 
             rep(c(0, 0, 1, 0), each = 25), rep(c(0, 0, 0, 1), each = 25), rnorm(n),
             rep(0:1, 50))
  B <- rbind(rnorm(J), seq(1, 5, length.out = J),
            rnorm(J), rnorm(J), rnorm(J), rnorm(J),
            rnorm(J), rnorm(J))
  for (k in 1:ncol(X)) {
    B[k, ] <- B[k, ] - radEmu:::pseudohuber_center(B[k, ], 0.1)
  }
  Y <- radEmu:::simulate_data(n = n,
                              J = J,
                              X = X,
                              B = B,
                              distn = "Poisson",
                              mean_z = 10)
  
  start <- proc.time()
  pl_fit_one <- emuFit_micro(X,
                             Y,
                             B = NULL,
                             maxit = 10000,
                             tolerance = 1e-6,
                             verbose= FALSE,
                             optimize_rows = FALSE)
  end <- proc.time() - start 
  start_refit <- proc.time() 
  pl_fit_two <- emuFit_micro(X,
                             Y,
                             B = NULL,
                             maxit = 10000,
                             tolerance =1e-6,
                             verbose= FALSE,
                             optimize_rows = TRUE)
  end_refit <- proc.time() - start_refit 
  expect_true(end_refit[3] < end[3])
  max_est_error1 <- max(pl_fit_one - B)
  max_est_error2 <- max(pl_fit_two - B)
  max_est_diff <- max(pl_fit_one - pl_fit_two)
  expect_true(max_est_error1 > max_est_error2)
  mean_est_error1 <- mean(sqrt((pl_fit_one - B)^2))
  mean_est_error2 <- mean(sqrt((pl_fit_two - B)^2))
  mean_est_diff <- mean(sqrt((pl_fit_one - pl_fit_two)^2))
  expect_true(mean_est_error1 > mean_est_error2)
  expect_true(mean_est_diff < 1e-3)
    
})
#
#
#

#
# test_that("ML fit to multiple regressors, large n, and excess-Poisson variance gives reasonable output", {
#   set.seed(4323)
#   n <- 800
#   X <- cbind(1,rep(c(0,1),each = n/2),rnorm(n),rnorm(n),rnorm(n))
#   z <- rnorm(n) +5
#   J <- 10
#   p <- 2
#
#   b0 <- rnorm(J)
#   b1 <- seq(1,10,length.out = J)
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = J, nrow = n)
#
#   for(i in 1:n){
#     for(j in 1:J){
#       temp_mean <- exp(X[i,1:2,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rnbinom(1,mu= temp_mean,size = 0.25)#rpois(1, lambda = temp_mean)
#     }
#   }
#   ml_fit <- emuFit_micro(X,
#                          Y,
#                          constraint_fn = function(x) mean(x),
#                          maxit = 200,
#                          tolerance = 1e-2)
#
#   # plot(b1-mean(b1),ml_fit[2,])
#   # abline(a = 0,b = 1,lty =2,col = "red")
#
#   expect_true(max(abs(ml_fit[2,] - (b1 - mean(b1))))<.5)
# })
#
#
# test_that("ML fit to multiple regressors, moderate n, moderate J, and excess-Poisson variance gives reasonable output", {
#   set.seed(4323)
#   n <- 40
#   X <- cbind(1,rep(c(0,1),each = n/2),rnorm(n),rnorm(n),rnorm(n))
#   z <- rnorm(n)
#   J <- 100
#   p <- 2
#
#   b0 <- rnorm(J)
#   b1 <- seq(1,10,length.out = J)
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = J, nrow = n)
#
#   for(i in 1:n){
#     for(j in 1:J){
#       temp_mean <- exp(X[i,1:2,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rnbinom(1,mu= temp_mean,size = 0.25)#rpois(1, lambda = temp_mean)
#     }
#   }
#   ml_fit <- emuFit_micro(X,
#                          Y,
#                          constraint_fn = function(x) mean(x),
#                          maxit = 200,
#                          tolerance = 1e-1,
#                          max_step = 0.5)
#
#   # plot(b1-mean(b1),ml_fit[2,])
#   # abline(a = 0,b = 1, lty = 2, col = "red")
#
#   b1_hat <- ml_fit[2,]
#   expect_true(abs(lm(b1_hat~b1)$coef[2] - 1) < 0.2)
#   # abline(a = 0,b = 1,lty =2,col = "red")
#
# })
#
# # t(ml_fit) %>% as.data.frame() %>% mutate(j = 1:J) %>%
# #   pivot_longer(-j) %>%
# #   ggplot() +
# #   geom_point(aes(x = j,y = value, color = name)) +
# #   geom_line(aes(x = j, y  = value, group = name, color = name)) +
# #   theme_bw()
#
#
#
# test_that("ML fit to simple example give reasonable output with J > n", {
#   set.seed(4323)
#   X <- cbind(1,rep(c(0,1),each = 20))
#   z <- rnorm(40) +8
#   b0 <- rnorm(80)
#   b1 <- seq(-5,5,length.out = 80)
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 80, nrow = 40)
#
#   for(i in 1:40){
#     for(j in 1:80){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rpois(1, lambda = temp_mean)
#     }
#   }
#   ml_fit <- emuFit_micro(X,
#                          Y,
#                          B = NULL,
#                          constraint_fn = function(x) mean(x),
#                          maxit = 500,
#                          tolerance = 0.01)
#
#
#
#
#
#   expect_true(max(abs(ml_fit[2,] - b1))<.1)
#
#   # X_cup <- X_cup_from_X(X,10)
#   # B_cup <- B_cup_from_B(ml_fit)
#
#
#
# })
#
#


