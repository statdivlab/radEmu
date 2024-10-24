
test_that("ML fit to simple example give reasonable output", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  z <- rnorm(40) +5
  J <- 10
  p <- 2
  n <- 40
  b0 <- rnorm(J)
  b1 <- seq(1,10,length.out = J)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = 40)

  for(i in 1:40){
    for(j in 1:J){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }
  ml_fit <- emuFit_micro(X,
             Y,
             constraint_fn = function(x) mean(x),
             maxit = 200,
             tolerance = 1e-3,
             verbose = FALSE)

  # plot(b1-mean(b1),ml_fit[2,])
  # abline(a = 0,b = 1,lty =2,col = "red")

  expect_true(max(abs(ml_fit[2,] - (b1 - mean(b1))))<.1)
})

test_that("With or without 'working_constraint' we get same results", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  z <- rnorm(40) +5
  J <- 10
  p <- 2
  n <- 40
  b0 <- rnorm(J)
  b1 <- seq(1,10,length.out = J)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = 40)

  for(i in 1:40){
    for(j in 1:J){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }
  ml_fit <- emuFit_micro(X,
                         Y,
                         constraint_fn = function(x) mean(x),
                         maxit = 200,
                         tolerance = 1e-6,
                         verbose= FALSE)
  ml_fit_direct <- emuFit_micro(X,
                                Y,
                                constraint_fn = function(x) mean(x),
                                maxit = 200,
                                use_working_constraint = FALSE,
                                tolerance = 1e-6,
                                verbose = FALSE)


  expect_equal(ml_fit,ml_fit_direct,tolerance = 1e-4)
})



test_that("PL fit with categorical predictor matches analytical form of MPLE in this case,
          and does NOT match MLE", {
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
                                             constraint_fn = function(x) mean(x),
                                             maxit = 200,
                                             tolerance = 1e-8,
                                             verbose= FALSE)
            
            ml_fit <- emuFit_micro(X,
                                   Y,
                                   B = matrix(rnorm(20),nrow = 2),
                                   constraint_fn = function(x) mean(x),
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
                                             constraint_fn = function(x) mean(x),
                                             maxit = 200,
                                             tolerance = 1e-8,
                                             verbose= FALSE)
            
            ml_fit <- emuFit_micro(X,
                                   Y,
                                   B = matrix(rnorm(20),nrow = 2),
                                   constraint_fn = function(x) mean(x),
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
  z <- rnorm(40) +5
  J <- 10
  p <- 2
  n <- 40
  b0 <- rnorm(J)
  b1 <- seq(1,10,length.out = J)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = 40)

  for(i in 1:40){
    for(j in 1:J){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }


  ml_fit <- emuFit_micro(X,
                         Y,
                         constraint_fn = function(x) mean(x),
                         maxit = 200,
                         tolerance = 1e-6,
                         verbose = FALSE)
  ml_fit_direct <- emuFit_micro(X,
                                Y,
                                constraint_fn = function(x) mean(x),
                                maxit = 200,
                                warm_start = FALSE,
                                tolerance = 1e-6,
                                verbose = FALSE)


  expect_equal(ml_fit,ml_fit_direct, tolerance = 1e-6)
})



test_that("We get a fit if we don't specify constraint", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  z <- rnorm(40) +5
  J <- 10
  p <- 2
  n <- 40
  b0 <- rnorm(J)
  b1 <- seq(1,10,length.out = J)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = 40)

  for(i in 1:40){
    for(j in 1:J){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }
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
  z <- rnorm(n) +8
  J <- 1000
  b0 <- rnorm(J)
  b1 <- seq(-5,5,length.out = J)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)

  for(i in 1:n){
    for(j in 1:J){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }


  ml_fit <- emuFit_micro(X,
                         Y,
                         constraint_fn = function(x) mean(x),
                         maxit = 500,
                         tolerance = 1e-3,
                         verbose = FALSE)


  expect_true(max(abs(ml_fit[2,] - b1))<.5)



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


