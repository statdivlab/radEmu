test_that("emuFit takes formulas and actually fits a model", {
  

  set.seed(343234)
  n <- 100
  X <- cbind(1,rep(c(0,1),each = n/2))
  J <- 10
  z <- rnorm(n) +8
  b0 <- rnorm(10)
  b1 <- 1:10
  b1 <- b1 - mean(b1)
  b1[5] <- pseudohuber_center(b1[-5],0.1)
  b0 <- b0 - mean(b0)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = 10, nrow = n)

  k_constr <- 2
  j_constr <- 5
  p <- 2

  constraint_fn <- function(x){ pseudohuber_center(x,0.1)}

  ##### Arguments to fix:

  constraint_grad_fn <- function(x){dpseudohuber_center_dx(x,0.1)}

  constraint_hess_fn <- function(x,ind_1,ind_2){hess_pseudohuber_center(x,0.1,ind_1,ind_2)}
  b[2,4] <- constraint_fn(b[2,-4])

  X_cup <- X_cup_from_X(X,J)

  Y[] <- 0
  for(i in 1:n){
    while(sum(Y[i,])==0){
      for(j in 1:10){
        temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
        Y[i,j] <- rnbinom(1,mu = temp_mean, size = 3)*rbinom(1,1,0.6)
      }
    }
  }

  covariates <- data.frame(group = X[,2])

  fitted_model <-
    emuFit(Y = Y,
         X = X,
         formula = ~group,
         data = covariates,
         verbose = FALSE,
         B_null_tol = 1e-2,
         tolerance = 0.01,
         tau = 2,
         run_score_test = TRUE,
         return_wald_p = TRUE)


  expect_true(all(fitted_model$coef$wald_p>0 & fitted_model$coef$wald_p<1))
  expect_true(inherits(fitted_model$B,"matrix"))
  expect_true(inherits(fitted_model$Y_augmented,"matrix"))
  expect_true(cor(fitted_model$B[2,],b1)>0.95)

})


test_that("We can use emuFit to create an emuFit object and then
call emuFit again *without* refitting model and it will return same results", {
  

  set.seed(94043234)
  n <- 100
  X <- cbind(1,rep(c(0,1),each = n/2))
  J <- 10
  z <- rnorm(n) +8
  b0 <- rnorm(10)
  b1 <- 1:10
  b1 <- b1 - mean(b1)
  b1[5] <- pseudohuber_center(b1[-5],0.1)
  b0 <- b0 - mean(b0)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = 10, nrow = n)
  Y[] <- 0
  for(i in 1:n){
    while(sum(Y[i,])==0){
      for(j in 1:10){
        temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
        # Y[i,j] <- rpois(1, lambda = temp_mean)
        Y[i,j] <- rnbinom(1,mu = temp_mean, size = 3)*rbinom(1,1,0.6)
      }
    }
  }

  covariates <- data.frame(group = X[,2])

  fitted_model <-
    emuFit(Y = Y,
           formula = ~group,
           data = covariates,
           run_score_test = FALSE,
           return_wald_p = FALSE)

  second_model <-
    emuFit(Y = Y,
           formula = ~group,
           data = covariates,
           refit = FALSE,
           run_score_test = FALSE,
           fitted_model = fitted_model)

  expect_identical(fitted_model,second_model)


})


test_that("emuFit takes formulas and actually fits a model (with score tests)", {
  

  set.seed(9983334)
  n <- 10
  X <- cbind(1,rep(c(0,1),each = n/2))
  J <- 5
  z <- rnorm(n) +8
  b0 <- rnorm(J)
  b1 <- 1:J
  b1 <- b1 - mean(b1)
  b1[5] <- pseudohuber_center(b1[-5],0.1)
  b0 <- b0 - mean(b0)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)

  k_constr <- 2
  j_constr <- 5
  p <- 2

  constraint_fn <- function(x){ pseudohuber_center(x,0.1)}

  ##### Arguments to fix:

  constraint_grad_fn <- function(x){dpseudohuber_center_dx(x,0.1)}

  constraint_hess_fn <- function(x,ind_1,ind_2){hess_pseudohuber_center(x,0.1,ind_1,ind_2)}
  b[2,4] <- constraint_fn(b[2,-4])

  X_cup <- X_cup_from_X(X,J)

  Y[] <- 0
  for(i in 1:n){
    while(sum(Y[i,])==0){
      for(j in 1:J){
        temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
        # Y[i,j] <- rpois(1, lambda = temp_mean)
        Y[i,j] <- rnbinom(1,mu = temp_mean, size = 3)*rbinom(1,1,0.6)
      }
    }
  }

  covariates <- data.frame(group = X[,2])

  fitted_model <-
    suppressMessages(
      emuFit(Y = Y,
           X = X,
           formula = ~group,
           tau = 1.2,
           data = covariates,
           run_score_test = TRUE,
           return_wald_p = TRUE,
           verbose = FALSE)
      )


  expect_true(all(fitted_model$coef$wald_p>0 & fitted_model$coef$wald_p<1))
  expect_true(inherits(fitted_model$B,"matrix"))
  expect_true(inherits(fitted_model$Y_augmented,"matrix"))
  expect_true(all(fitted_model$coef$pval>0 & fitted_model$coef$pval<1))


})

test_that("emuFit takes formulas and actually fits a model (with score tests using full model info)", {
  

  set.seed(9983334)
  n <- 10
  X <- cbind(1,rep(c(0,1),each = n/2))
  J <- 5
  z <- rnorm(n) +8
  b0 <- rnorm(J)
  b1 <- 1:J
  b1 <- b1 - mean(b1)
  b1[5] <- pseudohuber_center(b1[-5],0.1)
  b0 <- b0 - mean(b0)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)

  k_constr <- 2
  j_constr <- 5
  p <- 2

  constraint_fn <- function(x){ pseudohuber_center(x,0.1)}

  ##### Arguments to fix:

  constraint_grad_fn <- function(x){dpseudohuber_center_dx(x,0.1)}

  constraint_hess_fn <- function(x,ind_1,ind_2){hess_pseudohuber_center(x,0.1,ind_1,ind_2)}
  b[2,4] <- constraint_fn(b[2,-4])

  X_cup <- X_cup_from_X(X,J)

  Y[] <- 0
  for(i in 1:n){
    while(sum(Y[i,])==0){
      for(j in 1:J){
        temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
        # Y[i,j] <- rpois(1, lambda = temp_mean)
        Y[i,j] <- rnbinom(1,mu = temp_mean, size = 3)*rbinom(1,1,0.6)
      }
    }
  }

  covariates <- data.frame(group = X[,2])

  fitted_model <-
    suppressMessages(
      emuFit(Y = Y,
             X = X,
             formula = ~group,
             tau = 2,
             B_null_tol = 0.01,
             tolerance = 0.01,
             data = covariates,
             run_score_test = TRUE,
             return_wald_p = TRUE,
             use_fullmodel_info = TRUE,
             verbose = FALSE)
    )


  expect_true(all(fitted_model$coef$wald_p>0 & fitted_model$coef$wald_p<1))
  expect_true(inherits(fitted_model$B,"matrix"))
  expect_true(inherits(fitted_model$Y_augmented,"matrix"))
  expect_true(all(fitted_model$coef$pval>0 & fitted_model$coef$pval<1))


})


test_that("if return_both_score_pvals = TRUE, emuFit runs and returns two non-identical p-values", {
  

  set.seed(9983334)
  n <- 10
  X <- cbind(1,rep(c(0,1),each = n/2))
  J <- 5
  z <- rnorm(n) +8
  b0 <- rnorm(J)
  b1 <- 1:J
  b1 <- b1 - mean(b1)
  b1[5] <- pseudohuber_center(b1[-5],0.1)
  b0 <- b0 - mean(b0)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)

  k_constr <- 2
  j_constr <- 5
  p <- 2

  constraint_fn <- function(x){ pseudohuber_center(x,0.1)}

  ##### Arguments to fix:

  constraint_grad_fn <- function(x){dpseudohuber_center_dx(x,0.1)}

  constraint_hess_fn <- function(x,ind_1,ind_2){hess_pseudohuber_center(x,0.1,ind_1,ind_2)}
  b[2,4] <- constraint_fn(b[2,-4])

  X_cup <- X_cup_from_X(X,J)

  Y[] <- 0
  for(i in 1:n){
    while(sum(Y[i,])==0){
      for(j in 1:J){
        temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
        # Y[i,j] <- rpois(1, lambda = temp_mean)
        Y[i,j] <- rnbinom(1,mu = temp_mean, size = 3)*rbinom(1,1,0.6)
      }
    }
  }

  covariates <- data.frame(group = X[,2])

  fitted_model <-
    suppressMessages(
      emuFit(Y = Y,
             X = X,
             formula = ~group,
             tau = 1.2,
             data = covariates,
             run_score_test = TRUE,
             return_wald_p = TRUE,
             use_fullmodel_info = TRUE,
             verbose = FALSE,
             return_both_score_pvals = TRUE)
    )


  expect_true(is.numeric(fitted_model$coef$score_pval_full_info))
  expect_true(is.numeric(fitted_model$coef$score_pval_null_info))
  max_abs_logp_diff <- max(abs(log(fitted_model$coef$score_pval_null_info/
                                     fitted_model$coef$score_pval_full_info)))
  expect_true(max_abs_logp_diff>0)
  expect_true(max_abs_logp_diff<log(1.5))




})


test_that("if return_both_score_pvals = TRUE, we
          get same p-values as if we separately run emuFit with use_fullmodel_info = TRUE and = FALSE", {
  

  set.seed(334)
  n <- 10
  X <- cbind(1,rep(c(0,1),each = n/2))
  J <- 5
  z <- rnorm(n) +8
  b0 <- rnorm(J)
  b1 <- 1:J
  b1 <- b1 - mean(b1)
  b1[5] <- pseudohuber_center(b1[-5],0.1)
  b0 <- b0 - mean(b0)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)

  k_constr <- 2
  j_constr <- 5
  p <- 2

  constraint_fn <- function(x){ pseudohuber_center(x,0.1)}

  ##### Arguments to fix:

  constraint_grad_fn <- function(x){dpseudohuber_center_dx(x,0.1)}

  constraint_hess_fn <- function(x,ind_1,ind_2){hess_pseudohuber_center(x,0.1,ind_1,ind_2)}
  b[2,4] <- constraint_fn(b[2,-4])

  X_cup <- X_cup_from_X(X,J)

  Y[] <- 0
  for(i in 1:n){
    while(sum(Y[i,])==0){
      for(j in 1:J){
        temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
        # Y[i,j] <- rpois(1, lambda = temp_mean)
        Y[i,j] <- rnbinom(1,mu = temp_mean, size = 3)*rbinom(1,1,0.6)
      }
    }
  }

  covariates <- data.frame(group = X[,2])

  fitted_both <-
    suppressMessages(
      emuFit(Y = Y,
             X = X,
             formula = ~group,
             tau = 1.2,
             data = covariates,
             run_score_test = TRUE,
             return_wald_p = TRUE,
             use_fullmodel_info = TRUE,
             verbose = FALSE,
             return_both_score_pvals = TRUE)
    )

  fitted_null <-
    suppressMessages(
      emuFit(Y = Y,
             X = X,
             formula = ~group,
             tau = 1.2,
             data = covariates,
             run_score_test = TRUE,
             return_wald_p = TRUE,
             use_fullmodel_info = FALSE,
             verbose = FALSE,
             return_both_score_pvals = FALSE)
    )

  fitted_full <-
    suppressMessages(
      emuFit(Y = Y,
             X = X,
             formula = ~group,
             tau = 1.2,
             data = covariates,
             run_score_test = TRUE,
             return_wald_p = TRUE,
             use_fullmodel_info = TRUE,
             verbose = FALSE,
             return_both_score_pvals = FALSE)
    )


  expect_true(max(abs(fitted_both$coef$score_pval_full_info - fitted_full$coef$pval))
              ==0)
  expect_true(max(abs(fitted_both$coef$score_pval_null_info - fitted_null$coef$pval))
              ==0)

})
test_that("Difference between using full and null model info is not extreme", {
  

  set.seed(9983334)
  n <- 10
  X <- cbind(1,rep(c(0,1),each = n/2))
  J <- 5
  z <- rnorm(n) +8
  b0 <- rnorm(J)
  b1 <- 1:J
  b1 <- b1 - mean(b1)
  b1[5] <- pseudohuber_center(b1[-5],0.1)
  b0 <- b0 - mean(b0)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)

  k_constr <- 2
  j_constr <- 5
  p <- 2

  constraint_fn <- function(x){ pseudohuber_center(x,0.1)}

  ##### Arguments to fix:

  constraint_grad_fn <- function(x){dpseudohuber_center_dx(x,0.1)}

  constraint_hess_fn <- function(x,ind_1,ind_2){hess_pseudohuber_center(x,0.1,ind_1,ind_2)}
  b[2,4] <- constraint_fn(b[2,-4])

  X_cup <- X_cup_from_X(X,J)

  Y[] <- 0
  for(i in 1:n){
    while(sum(Y[i,])==0){
      for(j in 1:J){
        temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
        # Y[i,j] <- rpois(1, lambda = temp_mean)
        Y[i,j] <- rnbinom(1,mu = temp_mean, size = 3)*rbinom(1,1,0.6)
      }
    }
  }

  covariates <- data.frame(group = X[,2])

  fitted_model_full <-
    suppressMessages(
      emuFit(Y = Y,
             X = X,
             formula = ~group,
             tau = 1.2,
             data = covariates,
             run_score_test = TRUE,
             return_wald_p = TRUE,
             use_fullmodel_info = TRUE,
             verbose = FALSE)
    )

  fitted_model_null <-
    suppressMessages(
      emuFit(Y = Y,
             X = X,
             formula = ~group,
             tau = 1.2,
             data = covariates,
             run_score_test = TRUE,
             return_wald_p = TRUE,
             use_fullmodel_info = FALSE,
             verbose = FALSE)
    )


  expect_true(max(abs(log(fitted_model_full$coef$pval/fitted_model_null$coef$pval)))<
                log(1.5))




})

test_that("emuFit takes formulas and actually fits a model (no score tests) when J is large", {
  

  set.seed(894334)
  n <- 100
  X <- cbind(1,rep(c(0,1),each = n/2))
  J <- 500
  z <- rnorm(n) +8
  b0 <- rnorm(J)
  b1 <- seq(1,10,length.out = J)
  b1 <- b1 - mean(b1)
  b1[5] <- pseudohuber_center(b1[-5],0.1)
  b0 <- b0 - mean(b0)
  b <- rbind(b0,b1)
  Y <- matrix(NA,ncol = J, nrow = n)

  k_constr <- 2
  j_constr <- 5
  p <- 2

  constraint_fn <- function(x){ pseudohuber_center(x,0.1)}

  ##### Arguments to fix:

  constraint_grad_fn <- function(x){dpseudohuber_center_dx(x,0.1)}

  constraint_hess_fn <- function(x,ind_1,ind_2){hess_pseudohuber_center(x,0.1,ind_1,ind_2)}
  b[2,4] <- constraint_fn(b[2,-4])

  X_cup <- X_cup_from_X(X,J)

  Y[] <- 0
  for(i in 1:n){
    while(sum(Y[i,])==0){
      for(j in 1:J){
        temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
        # Y[i,j] <- rpois(1, lambda = temp_mean)
        Y[i,j] <- rnbinom(1,mu = temp_mean, size = 3)*rbinom(1,1,0.6)
      }
    }
  }

  covariates <- data.frame(group = X[,2])

  fitted_model <-
    emuFit(Y = Y,
           X = X,
           formula = ~group,
           data = covariates,
           tolerance = 0.01,
           run_score_test = FALSE,
           return_wald_p = TRUE)


  expect_true(all(fitted_model$coef$wald_p >= 0 & fitted_model$coef$wald_p <= 1))
  expect_true(inherits(fitted_model$B,"matrix"))
  expect_true(inherits(fitted_model$Y_augmented,"matrix"))
  expect_true(cor(fitted_model$B[2,],b1)>0.95)

})


# test_that("We can use emuFit to create an emuFit object and then
# call emuFit again to run score tests.", {
#   
#
#   set.seed(1243234)
#   n <- 100
#   X <- cbind(1,rep(c(0,1),each = n/2))
#   J <- 10
#   z <- rnorm(n) +8
#   b0 <- rnorm(10)
#   b1 <- 1:10
#   b1 <- b1 - mean(b1)
#   b1[5] <- pseudohuber_center(b1[-5],0.1)
#   b0 <- b0 - mean(b0)
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = 10, nrow = n)
#   Y[] <- 0
#   for(i in 1:n){
#     while(sum(Y[i,])==0){
#       for(j in 1:10){
#         temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#         # Y[i,j] <- rpois(1, lambda = temp_mean)
#         Y[i,j] <- rnbinom(1,mu = temp_mean, size = 3)*rbinom(1,1,0.6)
#       }
#     }
#   }
#
#   covariates <- data.frame(group = X[,2])
#
#   fitted_model <-
#     emuFit(Y = Y,
#            formula = ~group,
#            data = covariates,
#            run_score_test = FALSE,
#            return_wald_p = FALSE)
#
#   second_model <-
#     emuFit(Y = Y,
#            formula = ~group,
#            data = covariates,
#            refit = FALSE,
#            run_score_test = FALSE,
#            inner_maxit = 5,
#            verbose = TRUE,
#            tau = 1.1,
#            fitted_model = fitted_model)
#
#   expect_identical(fitted_model,second_model)
#
#
# })
