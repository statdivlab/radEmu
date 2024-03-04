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

# Y <- structure(c(534337, 0, 0, 0, 376, 41, 19, 103, 0, 0, 85, 0, 42794, 
#                  0, 0, 0, 95, 0, 0, 15, 0, 0, 0, 26, 0, 149, 0, 0, 0, 0, 0, 211, 
#                  0, 0, 0, 0, 0, 103, 0, 0, 0, 1372, 83, 337, 0, 0, 0, 0, 0, 53, 
#                  0, 0, 0, 0, 259, 0, 0, 0, 14, 0, 0, 0, 0, 193, 0, 0, 0, 0, 0, 
#                  0, 402, 0), dim = c(12L, 6L))
# X <- structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -2.49421928123597, 
#                  -0.579053917775617, -0.974555155010523, 0.237710056670222, 0.41240637454179, 
#                  1.12994912631468, 0.706485861932659, 0.588878125500377, 0.0834756145662259, 
#                  1.99483775157368, 0.227951737778031, -1.03963299361785), dim = c(12L, 
#                                                                                   2L))

# 
# Y <- structure(c(1748, 8286, 4096, 4289, 1122, 30007, 5087, 3841, 
#                  3059, 3105, 80, 32, 0, 20, 13, 0, 0, 30, 41, 54, 0, 124134, 43569, 
#                  122134, 15785, 99540, 0, 41104, 0, 0, 0, 0, 0, 572, 0, 1497, 
#                  0, 0, 314, 0, 0, 0, 0, 416, 920, 0, 0, 1931, 1279, 0, 0, 0, 2, 
#                  21, 49, 41, 0, 89, 0, 85, 1287, 1716, 0, 0, 1354, 8783, 3040, 
#                  6271, 2274, 0, 26431, 5186, 4147, 0, 0, 6450, 0, 0, 1483, 0, 
#                  0, 0, 0, 5936, 0, 0, 0, 33557, 11459, 0, 0, 4065, 0, 5391, 6721, 
#                  8997, 9225, 13951, 4061, 3871, 0, 0, 0, 0, 0, 0, 3954, 1338, 
#                  886, 426, 0, 0, 0, 496, 0, 709, 508, 840, 680, 0, 0, 4529, 2885, 
#                  0, 0, 13382, 11802, 0, 2144, 2622, 92214, 18326, 6183, 11737, 
#                  0, 0, 12808, 7604, 4684, 9348, 0, 3564, 2250, 0, 0, 20486, 0, 
#                  3442, 5133, 5103, 299927, 34273, 17407, 0, 17896, 149402, 42592, 
#                  0, 0, 25395, 1875, 21975, 1685, 0, 1654, 28670, 9331, 0, 4765, 
#                  0, 0, 0, 76158, 0, 85495, 247396, 51659, 93587, 0, 84277, 0, 
#                  0, 217, 0, 0, 0, 0, 0, 0, 0, 0, 76950, 0, 19911, 33042, 43690, 
#                  68014, 43885, 6086, 0), dim = c(20L, 10L))
# X <- structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
#                  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 
#                  1, 1), dim = c(20L, 2L))
covariates <- data.frame(group = X[,2])
# b <- structure(c(0.0231122313161752, -4.5, 2.34881909126062, -3.5, 
#                  -0.949623561962149, -2.5, -0.23942645176718, 0.556978794649417, 
#                  0.681191947270542, 0.499524350059682, -1.14509952873781, 0.5, 
#                  0.258373890958553, 1.5, 0.163812326430595, 2.5, 0.491510413255902, 
#                  3.5, -1.63267035802524, 4.5), dim = c(2L, 10L), dimnames = list(
#                    c("b0", "b1"), NULL))

test_that("emuFit takes formulas and actually fits a model", {
  
  ## test verbose = FALSE works
  expect_silent({
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
                           run_score_tests = TRUE, 
                           use_fullmodel_info = FALSE,
                           use_fullmodel_cov = FALSE,
                           return_both_score_pvals = FALSE)
  })
  
  
  
  ## emuFit takes formulas and actually fits a model (with score tests)
  
  expect_true(all(fitted_model$coef$wald_p>0 & fitted_model$coef$wald_p<1))
  expect_true(all(fitted_model$coef$pval>0 & fitted_model$coef$pval<1))
  expect_true(inherits(fitted_model$B,"matrix"))
  expect_true(inherits(fitted_model$Y_augmented,"matrix"))
  expect_true(cor(fitted_model$B[2,],b[2, ]) > 0.85) ## good enough for small sample size
  
  ## We can use emuFit to create an emuFit object and then
  ## call emuFit again *without* refitting model and it will return same results"
  second_model <- emuFit(Y = Y,
                         formula = ~group,
                         data = covariates,
                         refit = FALSE,
                         run_score_test = FALSE,
                         fitted_model = fitted_model)
  
  expect_identical(fitted_model$coef$estimate, second_model$coef$estimate)
  
  ##  emuFit takes formulas and actually fits a model (with score tests using full model info)
  
  fitted_model_use_fullmodel_info <- emuFit(Y = Y,
                                            X = X,
                                            formula = ~group,
                                            tau = 2,
                                            B_null_tol = 0.01,
                                            tolerance = 0.01,
                                            data = covariates,
                                            run_score_test = TRUE,
                                            return_wald_p = TRUE, ### diff
                                            use_fullmodel_info = TRUE, ### diff
                                            verbose = FALSE)
  
  
  expect_true(all(fitted_model_use_fullmodel_info$coef$wald_p>0 & fitted_model_use_fullmodel_info$coef$wald_p<1))
  expect_true(inherits(fitted_model_use_fullmodel_info$B,"matrix"))
  expect_true(inherits(fitted_model_use_fullmodel_info$Y_augmented,"matrix"))
  expect_true(all(fitted_model_use_fullmodel_info$coef$pval>0 & fitted_model_use_fullmodel_info$coef$pval<1))
  
  ## if return_both_score_pvals = TRUE, emuFit runs and returns two non-identical p-values
  
  ### TODO test this in a smaller example
  
  fitted_model_both <-  emuFit(Y = Y,
                               X = X,
                               formula = ~group,
                               tau = 1.2,
                               data = covariates,
                               run_score_test = TRUE,
                               return_wald_p = TRUE,
                               use_fullmodel_info = TRUE,
                               verbose = FALSE,
                               return_both_score_pvals = TRUE)
  
  ps_full <- fitted_model_both$coef$score_pval_full_info
  ps_null <- fitted_model_both$coef$score_pval_null_info
  expect_true(is.numeric(ps_full))
  expect_true(is.numeric(ps_null))
  expect_true(cor(ps_null, ps_full) > 0.95)
  expect_true(cor(ps_null, ps_full) < 1)
   
  # ## if return_both_score_pvals = TRUE, we
  # ## get same p-values as if we separately run emuFit with use_fullmodel_info = TRUE and = FALSE"
  # 
  # # fitted_model had use_fullmodel_info = FALSE
  # 
  # fitted_full <- emuFit(Y = Y,
  #                       X = X,
  #                       formula = ~group,
  #                       tau = 1.2,
  #                       data = covariates,
  #                       run_score_test = TRUE,
  #                       return_wald_p = TRUE,
  #                       use_fullmodel_info = TRUE,
  #                       verbose = FALSE,
  #                       return_both_score_pvals = FALSE)
  # 
  # expect_true(max(abs(fitted_model_both$coef$score_pval_full_info - fitted_full$coef$pval))==0)
  # expect_true(max(abs(fitted_model_both$coef$score_pval_null_info - fitted_model$coef$pval))==0)
  # 
  # ## Difference between using full and null model info is not extreme
  # expect_true(cor(fitted_full$coef$pval, fitted_model$coef$pval) < 1)
  # expect_true(cor(fitted_full$coef$pval, fitted_model$coef$pval) > 0.95)
  
})

test_that("emuFit takes cluster argument without breaking ",{
          expect_silent({
            fitted_model_cluster <- emuFit(Y = Y,
                                   X = X,
                                   formula = ~group,
                                   data = covariates,
                                   verbose = FALSE,
                                   B_null_tol = 1e-2,
                                   tolerance = 0.01,
                                   tau = 2,
                                   return_wald_p = FALSE,
                                   compute_cis = TRUE,
                                   run_score_tests = TRUE, 
                                   use_fullmodel_info = FALSE,
                                   use_fullmodel_cov = FALSE,
                                   return_both_score_pvals = FALSE,
                                   cluster = rep(1:3,each = 4))
          })
  
  expect_silent({
    fitted_model_nocluster <- emuFit(Y = Y,
                                   X = X,
                                   formula = ~group,
                                   data = covariates,
                                   verbose = FALSE,
                                   B_null_tol = 1e-2,
                                   tolerance = 0.01,
                                   tau = 2,
                                   return_wald_p = FALSE,
                                   compute_cis = TRUE,
                                   run_score_tests = TRUE, 
                                   use_fullmodel_info = FALSE,
                                   use_fullmodel_cov = FALSE,
                                   return_both_score_pvals = FALSE)
  })
  
  expect_true(all(fitted_model_nocluster$coef$estimate == fitted_model_cluster$coef$estimate))
  
  expect_true(all(fitted_model_nocluster$coef$se != fitted_model_cluster$coef$se))
  
  expect_true(all(fitted_model_cluster$coef$se != fitted_model_nocluster$coef$se))
  
          
          
          
          ## emuFit takes formulas and actually fits a model (with score tests)
          
          expect_true(all(fitted_model_cluster$coef$wald_p>0 & fitted_model_cluster$coef$wald_p<1))
          expect_true(all(fitted_model_cluster$coef$pval>0 & fitted_model_cluster$coef$pval<1))
          expect_true(inherits(fitted_model_cluster$B,"matrix"))
          expect_true(inherits(fitted_model_cluster$Y_augmented,"matrix"))
          expect_true(cor(fitted_model_cluster$B[2,],b[2, ]) > 0.85)
          
})

b0 <- rnorm(J)
b1 <- seq(1,5,length.out = J)
b1 <- b1 - mean(b1)
b1[3:4] <- 0
b <- rbind(b0,b1)

test_that("GEE with cluster covariance gives plausible type 1 error ",{
  skip("Skipping -- test requires fitting models to 100 simulated datasets.")
  set.seed(44022)
  nsim <- 100
  cluster <- rep(1:4, each = 3)
  results <- data.frame(sim = rep(1:nsim,each = 2),
                        category_num = rep(3:4,nsim),
                        estimate = numeric(2*nsim),
                        score_stat = numeric(2*nsim),
                        pval = numeric(2*nsim))[-(1:(2*nsim)),]
  results_noGEE <- results
  for(sim in 1:nsim){
    print(sim)
    X <- cbind(1,rnorm(n))
    covariates <- data.frame(group = X[,2])
    Y <- matrix(NA,ncol = J, nrow = n)
  
  cluster_effs <- lapply(1:4,
                         function(i)
                           log(matrix(rexp(2*J),nrow= 2)))
  
  for(i in 1:n){
    Y[i,] <- 0
    while(sum(Y[i,])==0){
    for(j in 1:J){
      temp_mean <- exp(X[i,,drop = FALSE]%*%(b[,j,drop = FALSE] + 
                                               cluster_effs[[ cluster[i] ]][,j]) + z[i])
      Y[i,j] <- rnbinom(1, mu= temp_mean,size = 5)*rbinom(1,1,0.8)
    }}
  }
  
  # expect_silent({
    fitted_model_cluster <- emuFit(Y = Y,
                                   X = X,
                                   formula = ~group,
                                   data = covariates,
                                   verbose = FALSE,
                                   B_null_tol = 1e-2,
                                   tolerance = 0.01,
                                   tau = 1.2,
                                   return_wald_p = FALSE,
                                   compute_cis = TRUE,
                                   run_score_tests = TRUE, 
                                   use_fullmodel_info = FALSE,
                                   use_fullmodel_cov = FALSE,
                                   return_both_score_pvals = FALSE,
                                   test_kj = data.frame(k = c(2,2),
                                                        j = c(3,4)),
                                   cluster = cluster)
    
    fitted_model_nocluster <- emuFit(Y = Y,
                                   X = X,
                                   formula = ~group,
                                   data = covariates,
                                   verbose = FALSE,
                                   B_null_tol = 1e-2,
                                   tolerance = 0.01,
                                   tau = 1.2,
                                   return_wald_p = FALSE,
                                   compute_cis = TRUE,
                                   run_score_tests = TRUE, 
                                   use_fullmodel_info = FALSE,
                                   use_fullmodel_cov = FALSE,
                                   return_both_score_pvals = FALSE,
                                   test_kj = data.frame(k = c(2,2),
                                                        j = c(3,4)))
  # })
  
  filtered_coef <- fitted_model_cluster$coef[!is.na(fitted_model_cluster$coef$pval),
                                             c("category_num",
                                               "estimate",
                                               "score_stat",
                                               "pval")]
  results <- rbind(results,
                   cbind(data.frame("sim" = rep(sim,2)),
                         filtered_coef))
  
  filtered_coef_noGEE <- fitted_model_nocluster$coef[!is.na(fitted_model_nocluster$coef$pval),
                                             c("category_num",
                                               "estimate",
                                               "score_stat",
                                               "pval")]
  results_noGEE <- rbind(results_noGEE,
                   cbind(data.frame("sim" = rep(sim,2)),
                         filtered_coef_noGEE))
  
  
  }
  
  #expect somewhat conservative inference for score test with correct clustering
  expect_gte(mean(results$pval<=0.05), 0.05)
  #expect fairly anti-conservative inference for score test without clustering
  expect_gte(mean(results_noGEE$pval<=0.05), 0.10)
  
})

# test_that("emuFit takes formulas and actually fits a model (no score tests) when J is large", {
#   
#   b1 <- 1:10
# b1 <- b1 - mean(b1)
# b1[5] <- pseudohuber_center(b1[-5],0.1)
# 
# Y <- simulate_data(n=10, J=10, b0=rnorm(10), distn="Poisson", b1=b1, mean_count_before_ZI=500)

#   set.seed(894334)
#   n <- 100
#   X <- cbind(1,rep(c(0,1),each = n/2))
#   J <- 500
#   z <- rnorm(n) +8
#   b0 <- rnorm(J)
#   b1 <- seq(1,10,length.out = J)
#   b1 <- b1 - mean(b1)
#   b1[5] <- pseudohuber_center(b1[-5],0.1)
#   b0 <- b0 - mean(b0)
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = J, nrow = n)
#   
#   k_constr <- 2
#   j_constr <- 5
#   p <- 2
#   
#   constraint_fn <- function(x){ pseudohuber_center(x,0.1)}
#   
#   ##### Arguments to fix:
#   
#   constraint_grad_fn <- function(x){dpseudohuber_center_dx(x,0.1)}
#   
#   constraint_hess_fn <- function(x,ind_1,ind_2){hess_pseudohuber_center(x,0.1,ind_1,ind_2)}
#   b[2,4] <- constraint_fn(b[2,-4])
#   
#   X_cup <- X_cup_from_X(X,J)
#   
#   Y[] <- 0
#   for(i in 1:n){
#     while(sum(Y[i,])==0){
#       for(j in 1:J){
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
#            X = X,
#            formula = ~group,
#            data = covariates,
#            tolerance = 0.01,
#            run_score_test = FALSE,
#            return_wald_p = TRUE)
#   
#   
#   expect_true(all(fitted_model$coef$wald_p >= 0 & fitted_model$coef$wald_p <= 1))
#   expect_true(inherits(fitted_model$B,"matrix"))
#   expect_true(inherits(fitted_model$Y_augmented,"matrix"))
#   expect_true(cor(fitted_model$B[2,],b1)>0.95)
#   
# })
