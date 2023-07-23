test_that("Test score test with Poisson data under alternative", {
  set.seed(4323)
  X <- cbind(1,rep(c(0,1),each = 20))
  z <- rnorm(40) +8
  b0 <- rnorm(10)
  b1 <- 1:10
  b1 <- b1 - mean(b1)
  b0 <- b0 - mean(b0)
  b <- rbind(b0,b1)

  pvals <- numeric(100)
  for(iter in 1:100){
    print(iter)
  Y <- matrix(NA,ncol = 10, nrow = 40)

  for(i in 1:40){
    for(j in 1:10){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      Y[i,j] <- rpois(1, lambda = temp_mean)
    }
  }

  pvals[iter] <- micro_score_test(Y= Y,
                   X = X,
                   B= NULL,
                   constraint_fn = mean,
                   constraint_type = "mean",
                   null_k = 2,
                   null_j = 4)$pval


  }

  hist(pvals)

  expect_true(max(pvals)<2e-3)

})

test_that("Test score test with Poisson data under null, J > n", {
  set.seed(4323)
  n <- 10
  J <- 20
  X <- cbind(1,rep(c(-.5,.5),each = n/2))
  z <- rnorm(n) +3
  b0 <- rnorm(J)
  b1 <- seq(-1.5,1.5,length.out = J)
  b1[J/2 + 0:1] <- 0
  b <- rbind(b0,b1)

  nsim = 1e4
  pvals <- numeric(nsim)
  score_stats <- numeric(nsim)
  for(iter in 1:nsim){
    print(iter)
    Y <- matrix(NA,ncol = J, nrow = n)

    for(i in 1:n){
      for(j in 1:J){
        temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
        # Y[i,j] <- rpois(1, lambda = temp_mean)
        Y[i,j] <- rnbinom(1, mu = temp_mean,size = 1.5)
      }
    }

    score_result <- try(micro_score_test(Y= Y,
                                     X = X,
                                     B= NULL,
                                     constraint_fn = mean,
                                     null_k = 2,
                                     null_j = J/2,
                                     tolerance = 1e-5,
                                     constraint_type = "mean"))


    pvals[iter] <- score_result$pval
    # print(pvals[iter])

    score_stats[iter] <- score_result$score_stat
    if(iter %%10 ==0){
    # hist(score_stats[1:iter],xlim = c(-10,10),freq= FALSE,breaks= 50)
    # x = seq(-10,10,by = .01)
    # lines(x,dt(x,n - 1),lty= 2)
      hist(score_stats[1:iter],xlim = c(0,10),freq= FALSE,breaks= 50)
    x <- seq(0,10,by = 0.01)
    lines(x,dchisq(x,1))
    print(signif(mean(pvals[1:iter]<=0.05)))
}

  }

hist(pvals,freq = FALSE,breaks = 100)
 mean(pvals <=0.05)

hist(score_stats,breaks =200,freq=FALSE)
x = seq(-10,10,by = .01)
lines(x,dt(x,n - 1),lty= 2,col="red")

qqplot(score_stats,rt(1e5,n-1))
abline(a= 0, b=1, col  ="red")
sd(score_stats)


})

test_that("Test score test with Poisson data under null, J > n, huber constraint", {
  set.seed(4323)
  n <- 10
  J <- 20
  X <- cbind(1,rep(c(-.5,.5),each = n/2))
  z <- rnorm(n) +3
  b0 <- rnorm(J)
  b1 <- seq(-1.5,1.5,length.out = J)
  b1[J/2 + 0:1] <- 0
  b <- rbind(b0,b1)

  nsim = 1e4
  pvals <- numeric(nsim)
  score_stats <- numeric(nsim)
  for(iter in 1:nsim){
    print(iter)
    Y <- matrix(NA,ncol = J, nrow = n)

    for(i in 1:n){
      for(j in 1:J){
        temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
        # Y[i,j] <- rpois(1, lambda = temp_mean)
        Y[i,j] <- rnbinom(1, mu = temp_mean,size = 1.5)
      }
    }

    score_result <- try(micro_score_test(Y= Y,
                                         X = X,
                                         B= NULL,
                                         constraint_fn = huber_center,
                                         null_k = 2,
                                         null_j = J/2,
                                         tolerance = 1e-4,
                                         constraint_type = "huber",
                                         huber_param = 1))


    pvals[iter] <- score_result$pval
    # print(pvals[iter])

    score_stats[iter] <- score_result$score_stat
    if(iter %%10 ==0){
      # hist(score_stats[1:iter],xlim = c(-10,10),freq= FALSE,breaks= 50)
      # x = seq(-10,10,by = .01)
      # lines(x,dt(x,n - 1),lty= 2)
      hist(score_stats[1:iter],xlim = c(0,10),freq= FALSE,breaks= 50)
      x <- seq(0,10,by = 0.01)
      lines(x,dchisq(x,1))
      print(signif(mean(pvals[1:iter]<=0.05)))
    }

  }

  hist(pvals,freq = FALSE,breaks = 100)
  mean(pvals <=0.05)

  hist(score_stats,breaks =200,freq=FALSE)
  x = seq(-10,10,by = .01)
  lines(x,dt(x,n - 1),lty= 2,col="red")

  qqplot(score_stats,rt(1e5,n-1))
  abline(a= 0, b=1, col  ="red")
  sd(score_stats)


})

test_that("Test score test with Poisson data under null, larger n", {
  set.seed(4323)
  n <- 20
  J <- 80
  X <- cbind(1,rep(c(-.5,.5),each = n/2))
  z <- rnorm(n) +3
  b0 <- rnorm(J)
  b1 <- seq(-1.5,1.5,length.out = J)
  b1[J/2 + 0:1] <- 0
  b <- rbind(b0,b1)

  nsim = 1e4
  pvals <- numeric(nsim)
  score_stats <- numeric(nsim)
  for(iter in 1:nsim){
    print(iter)
    Y <- matrix(NA,ncol = J, nrow = n)

    for(i in 1:n){
      for(j in 1:J){
        temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
        Y[i,j] <- rpois(1, lambda = temp_mean)
        # Y[i,j] <- rnbinom(1, mu = temp_mean,size = 1.5)
      }
    }

    score_result <- try(micro_score_test(Y= Y,
                                         X = X,
                                         B= NULL,
                                         constraint_fn = mean,
                                         null_k = 2,
                                         null_j = J/2,
                                         tolerance = 1e-5,
                                         type = "correct_score",
                                         cov_est_method = "plugin",
                                         t_se_method= "from_cov",
                                         loo = FALSE))


    pvals[iter] <- score_result$pval
    # print(pvals[iter])

    score_stats[iter] <- score_result$score_stat
    if(iter %%10 ==0){
      # hist(score_stats[1:iter],xlim = c(-10,10),freq= FALSE,breaks= 50)
      # x = seq(-10,10,by = .01)
      # lines(x,dt(x,n - 1),lty= 2)
      hist(score_stats[1:iter],xlim = c(0,10),freq= FALSE,breaks= 50)
      x <- seq(0,10,by = 0.01)
      lines(x,dchisq(x,1))
      print(signif(mean(pvals[1:iter]<=0.05)))
    }

  }

hist(pvals,freq = FALSE,breaks = 100)
  mean(pvals <=0.05)

  hist(score_stats,breaks =200,freq=FALSE)
  x = seq(-10,10,by = .01)
  lines(x,dt(x,n - 1),lty= 2,col="red")

  qqplot(score_stats,rt(1e5,n-1))
  abline(a= 0, b=1, col  ="red")
  sd(score_stats)


})

test_that("Test score test with Poisson data under null, J > n", {
  set.seed(4323)
  n <- 6
  J <- 20
  X <- cbind(1,rep(c(-.5,.5),each = n/2))
  z <- rnorm(n) +3
  b0 <- rnorm(J)
  b1 <- seq(-1.5,1.5,length.out = J)
  b1[J/2 + 0:1] <- 0
  b <- rbind(b0,b1)

  results <- vector(4,mode = "list")

  nsim = 1e4
  pvals <- numeric(nsim)
  score_stats <- numeric(nsim)
  for(iter in 1:nsim){
    print(iter)
    Y <- matrix(NA,ncol = J, nrow = n)

    for(i in 1:n){
      for(j in 1:J){
        temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
        # Y[i,j] <- rpois(1, lambda = temp_mean)
        Y[i,j] <- rnbinom(1, mu = temp_mean,size = 2)
      }
    }

    score_result <- try(micro_score_test(Y= Y,
                                         X = X,
                                         B= NULL,
                                         constraint_fn = mean,
                                         null_k = 2,
                                         null_j = J/2,
                                         tolerance = 1e-5,
                                         maxit = 50,
                                         type = "t",
                                         cov_est_method = "crossproduct",
                                         t_se_method= "from_cov",
                                         loo = FALSE))


    pvals[iter] <- score_result$pval
    # print(pvals[iter])

    score_stats[iter] <- score_result$score_stat
    if(iter %%10 ==0){
      # hist(score_stats[1:iter],xlim = c(-10,10),freq= FALSE,breaks= 50)
      # x = seq(-10,10,by = .01)
      # lines(x,dt(x,n - 1),lty= 2)
      hist(score_stats[1:iter],xlim = c(0,10),freq= FALSE,breaks= 50)
      x <- seq(0,10,by = 0.01)
      lines(x,dchisq(x,1),col = "red")
      print(signif(mean(pvals[1:iter]<=0.05)))
    }

  }

  results[[1]] <- list("score_stats" = score_stats,
                       "pvals" = pvals)

  pvals <- numeric(nsim)
  score_stats <- numeric(nsim)
  for(iter in 1:nsim){
    print(iter)
    Y <- matrix(NA,ncol = J, nrow = n)

    for(i in 1:n){
      for(j in 1:J){
        temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
        # Y[i,j] <- rpois(1, lambda = temp_mean)
        Y[i,j] <- rnbinom(1, mu = temp_mean,size = 1.5)
      }
    }

    score_result <- try(micro_score_test(Y= Y,
                                         X = X,
                                         B= NULL,
                                         constraint_fn = mean,
                                         null_k = 2,
                                         null_j = J/2,
                                         tolerance = 1e-5,
                                         type = "t",
                                         cov_est_method = "crossproduct",
                                         t_se_method= "from_cov",
                                         loo = FALSE))


    pvals[iter] <- score_result$pval
    # print(pvals[iter])

    score_stats[iter] <- score_result$score_stat
    if(iter %%10 ==0){
      # hist(score_stats[1:iter],xlim = c(-10,10),freq= FALSE,breaks= 50)
      # x = seq(-10,10,by = .01)
      # lines(x,dt(x,n - 1),lty= 2)
      hist(score_stats[1:iter],xlim = c(0,10),freq= FALSE,breaks= 50)
      x <- seq(0,10,by = 0.01)
      lines(x,dchisq(x,1))
      print(signif(mean(pvals[1:iter]<=0.05)))
    }

  }

  results[[2]] <- list("score_stats" = score_stats,
                       "pvals" = pvals)
  set.seed(4323)
  n <- 10
  J <- 80
  X <- cbind(1,rep(c(-.5,.5),each = n/2))
  z <- rnorm(n) +3
  b0 <- rnorm(J)
  b1 <- seq(-1.5,1.5,length.out = J)
  b1[J/2 + 0:1] <- 0
  b <- rbind(b0,b1)
  pvals <- numeric(nsim)
  score_stats <- numeric(nsim)
  for(iter in 1:nsim){
    print(iter)
    Y <- matrix(NA,ncol = J, nrow = n)

    for(i in 1:n){
      for(j in 1:J){
        temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
        Y[i,j] <- rpois(1, lambda = temp_mean)
        # Y[i,j] <- rnbinom(1, mu = temp_mean,size = 1.5)
      }
    }

    score_result <- try(micro_score_test(Y= Y,
                                         X = X,
                                         B= NULL,
                                         constraint_fn = mean,
                                         null_k = 2,
                                         null_j = J/2,
                                         tolerance = 1e-5,
                                         type = "t",
                                         cov_est_method = "crossproduct",
                                         t_se_method= "from_cov",
                                         loo = FALSE))


    pvals[iter] <- score_result$pval
    # print(pvals[iter])

    score_stats[iter] <- score_result$score_stat
    if(iter %%10 ==0){
      # hist(score_stats[1:iter],xlim = c(-10,10),freq= FALSE,breaks= 50)
      # x = seq(-10,10,by = .01)
      # lines(x,dt(x,n - 1),lty= 2)
      hist(score_stats[1:iter],xlim = c(0,10),freq= FALSE,breaks= 50)
      x <- seq(0,10,by = 0.01)
      lines(x,dchisq(x,1))
      print(signif(mean(pvals[1:iter]<=0.05)))
    }

  }

  results[[3]] <- list("score_stats" = score_stats,"pvals" = pvals)

  pvals <- numeric(nsim)
  score_stats <- numeric(nsim)
  for(iter in 1:nsim){
    print(iter)
    Y <- matrix(NA,ncol = J, nrow = n)

    for(i in 1:n){
      for(j in 1:J){
        temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
        # Y[i,j] <- rpois(1, lambda = temp_mean)
        Y[i,j] <- rnbinom(1, mu = temp_mean,size = 1.5)
      }
    }

    score_result <- try(micro_score_test(Y= Y,
                                         X = X,
                                         B= NULL,
                                         constraint_fn = mean,
                                         null_k = 2,
                                         null_j = J/2,
                                         tolerance = 1e-5,
                                         type = "t",
                                         cov_est_method = "crossproduct",
                                         t_se_method= "from_cov",
                                         loo = FALSE))


    pvals[iter] <- score_result$pval
    # print(pvals[iter])

    score_stats[iter] <- score_result$score_stat
    if(iter %%10 ==0){
      # hist(score_stats[1:iter],xlim = c(-10,10),freq= FALSE,breaks= 50)
      # x = seq(-10,10,by = .01)
      # lines(x,dt(x,n - 1),lty= 2)
      hist(score_stats[1:iter],xlim = c(0,10),freq= FALSE,breaks= 50)
      x <- seq(0,10,by = 0.01)
      lines(x,dchisq(x,1))
      print(signif(mean(pvals[1:iter]<=0.05)))
    }

  }

  results[[4]] <- list("score_stats" = score_stats,"pvals" = pvals)

})

data.frame("J" = rep(c(20,80),each= 2),
           "n" = 10,
           "Distribution" = rep(c("Poisson","NB with size = 1.5")),
           "empirical_alpha"= sapply(1:4,function(i) mean(results[[i]]$pval<=0.05))) %>%
  ggplot() +
  geom_point(aes(x = J,y = empirical_alpha,color = Distribution)) +
  ylim(c(0,.1)) +
  geom_abline(aes(intercept = 0.05,slope= 0),linetype = 2)+
  xlab("Number of Taxa") +
  ylab("Empirical Rejection Rate at 0.05 Level") +
  theme_bw()

rbind(
  data.frame("Distribution" = "Poisson",
             "pvals" = results[[1]]$pvals,
             "J" = 20),
  data.frame("Distribution" = "NB with size = 1.5",
             "pvals" = results[[2]]$pvals,
             "J" = 20),
  data.frame("Distribution" = "Poisson",
             "pvals" = results[[3]]$pvals,
             "J" = 80),
  data.frame("Distribution" = "NB with size = 1.5",
             "pvals" = results[[4]]$pvals,
             "J" = 80)
) %>%
  ggplot(aes(sample = pvals)) +
  stat_qq(distribution = qunif,size = 0.5,alpha = 0.01) +
  geom_abline(aes(intercept = 0, slope = 1),linetype = 2, color= "red") +
  facet_grid(J~Distribution) +
  theme_bw()


test_that("Test score test with Poisson data under null", {
  set.seed(4323)
  n <- 20
  J <- 10
  X <- cbind(1,rep(c(-.5,.5),each = n/2))
  z <- rnorm(n) +1
  b0 <- rnorm(J)
  b1 <- seq(-1.5,1.5,length.out = J)
  b1[J/2 + 0:1] <- 0
  b <- rbind(b0,b1)

  nsim = 1e4
  pvals <- numeric(nsim)
  score_stats <- numeric(nsim)
  for(iter in 1:nsim){
    print(iter)
    Y <- matrix(NA,ncol = J, nrow = n)

    for(i in 1:n){
      for(j in 1:J){
        temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
        Y[i,j] <- rpois(1, lambda = temp_mean)
        # Y[i,j] <- rnbinom(1, mu = temp_mean,size = 1.5)
      }
    }

    pf <- suppressMessages(emuFit_micro_penalized(X[rowSums(Y)>0,],
                  Y[rowSums(Y)>0,],
                  B = NULL,
                  constraint_fn = function(x) x[J],
                  maxit = 500,
                  tolerance = 1e-3,
                  collect_iterations = FALSE))

    score_result <- try(micro_score_test(Y= pf$Y_augmented,
                                         X = X[rowSums(Y)>0,],
                                         B= NULL,
                                         constraint_fn = mean,
                                         null_k = 2,
                                         null_j = J/2,
                                         tolerance = 1e-5,
                                         type = "correct_score",
                                         cov_est_method = "plugin",
                                         t_se_method= "from_cov",
                                         loo = FALSE))


    pvals[iter] <- score_result$pval
    # print(pvals[iter])

    score_stats[iter] <- score_result$score_stat
    if(iter %%10 ==0){
      # hist(score_stats[1:iter],xlim = c(-10,10),freq= FALSE,breaks= 50)
      # x = seq(-10,10,by = .01)
      # lines(x,dt(x,n - 1),lty= 2)
      hist(score_stats[1:iter],xlim = c(0,10),freq= FALSE,breaks= 50)
      x <- seq(0,10,by = 0.01)
      lines(x,dchisq(x,1))
      print(signif(mean(pvals[1:iter]<=0.05)))
    }

  }

    hist(pvals,freq = FALSE,breaks = 100)
  mean(pvals <=0.05)

  hist(score_stats,breaks =200,freq=FALSE)
  x = seq(-10,10,by = .01)
  lines(x,dt(x,n - 1),lty= 2,col="red")

  qqplot(score_stats,rt(1e5,n-1))
  abline(a= 0, b=1, col  ="red")
  sd(score_stats)


})


