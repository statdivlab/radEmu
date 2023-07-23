test_that("SEs seem reasonable in simple case", {
  set.seed(4323)
  n <- 20
  X <- cbind(1,rep(c(0,1),each = n/2))
  z <- rnorm(n) +4
  J <- 10
  p <- 2

  b0 <- rnorm(J)
  b1 <- 1:J
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)

  for(i in 1:n){
    for(j in 1:J){
      temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
      # Y[i,j] <- rpois(1, lambda = temp_mean)
      Y[i,j] <- rnbinom(1,size = 1.5,mu = temp_mean)
    }
  }
  ml_fit <- emuFit_micro(X,
                         Y,
                         B = NULL,
                         constraint_fn = function(x) mean(x),
                         maxit = 200,
                         tolerance = 1e-1)

  oi <- observed_info(X,Y,ml_fit)/n
  oi_ <- MASS::ginv(oi)

  V <- robust_V(X,Y,ml_fit)
  B_var_est <- oi_%*%V%*%oi_

  B_ses <- numeric(p*J)
  counter <- 1
  for(j in 1:J){
    for(k in 1:2){
      contrast_vec <- matrix(0,nrow = p*J,ncol = 1)
      contrast_vec[((1:(p*J)) - k )%%2 ==0,] <- -1/J
      contrast_vec[(j - 1)*p + k,] <- contrast_vec[(j - 1)*p + k,] + 1
      B_ses[counter] <- sqrt(as.numeric(t(contrast_vec)%*%B_var_est%*%contrast_vec))/sqrt(n)
      counter <- counter + 1
    }
  }

  long_truth <- as.numeric(as.matrix(B_cup_from_B(rbind(b0 - mean(b0),b1 - mean(b1)))))

  long_est <- as.numeric(as.matrix(B_cup_from_B(ml_fit)))

  data.frame("truth" = long_truth,
             "estimate" = long_est,
             "se" = B_ses,
             "k" = rep(1:2,J),
             "j" = rep(1:J, each = p)) %>%
    mutate(lower = estimate - qnorm(.975)*se,
           upper = estimate + qnorm(.975)*se) %>%
    ggplot() +
    geom_point(aes(x = interaction(k,j),y = estimate - truth)) +
    geom_errorbar(aes(x = interaction(k,j),ymin = lower - truth, ymax = upper - truth),width = 0) +
    theme_bw()

})
