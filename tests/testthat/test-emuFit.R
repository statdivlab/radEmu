set.seed(11)
J <- 6
n <- 12
X <- cbind(1,rnorm(n))
b0 <- rnorm(J)
b1 <- seq(1,5,length.out = J) -
  mean(seq(1,5,length.out = J))
b <- rbind(b0, b1)
Y <- radEmu:::simulate_data(n = n,
                            J = J,
                            X = X,
                            b0 = b0,
                            b1 = b1,
                            distn = "ZINB",
                            zinb_size = 2,
                            zinb_zero_prop = 0.2,
                            mean_z = 5)

#To ensure the messages about lack of row names do not show in the tests
rownames(X) <- paste0("Sample_",1:12)
rownames(Y) <- paste0("Sample_",1:12)

covariates <- data.frame(group = X[,2])

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
                           return_both_score_pvals = FALSE,
                           test_kj = data.frame(k = 2, j = 1:6))
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
                         run_score_tests= FALSE,
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
                                            run_score_tests= TRUE,
                                            return_wald_p = TRUE, ### diff
                                            use_fullmodel_info = TRUE, ### diff
                                            verbose = FALSE,
                                            test_kj = data.frame(k = 2, j = 1:6))
  
  
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
                               run_score_tests= TRUE,
                               return_wald_p = TRUE,
                               use_fullmodel_info = TRUE,
                               verbose = FALSE,
                               return_both_score_pvals = TRUE,
                               test_kj = data.frame(k = 2, j = 1:6))
  
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
  #                       run_score_tests= TRUE,
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
                                   cluster = rep(1:3,each = 4),
                                   test_kj = data.frame(k = 2, j = 1:6))
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
                                   return_both_score_pvals = FALSE,
                                   test_kj = data.frame(k = 2, j = 1:6))
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
    
    X <- cbind(1,rnorm(12))
    covariates <- data.frame(group = X[,2])
    
    b1 <- seq(1,5,length.out = 6)
    b1 <- b1 - mean(b1)
    b1[3:4] <- 0
    
    Y <- radEmu:::simulate_data(n = 12, J = 6,
                                X = X,
                                b0 = rnorm(6),
                                b1 = b1,
                                distn = "ZINB",
                                zinb_size = 2,
                                zinb_zero_prop = 0.2,
                                mean_z = 5,
                                cluster = cluster)
    
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
# Y <- radEmu:::simulate_data(n=10, J=10, b0=rnorm(10), distn="Poisson", b1=b1, mean_z=500)

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
#            run_score_tests= FALSE,
#            return_wald_p = TRUE)
#   
#   
#   expect_true(all(fitted_model$coef$wald_p >= 0 & fitted_model$coef$wald_p <= 1))
#   expect_true(inherits(fitted_model$B,"matrix"))
#   expect_true(inherits(fitted_model$Y_augmented,"matrix"))
#   expect_true(cor(fitted_model$B[2,],b1)>0.95)
#   
# })

test_that("emuFit runs without penalty", {
  
  expect_silent({
    fitted_model <- emuFit(Y = Y,
                           X = X,
                           penalize = FALSE,
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
                           test_kj = data.frame(k = 2, j = 1:6))
  })
})

test_that("emuFit runs with just intercept model", {
  
  expect_message({
    fitted_model <- emuFit(Y = Y,
                           formula = ~1,
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
                           test_kj = data.frame(k = 1, j = 1:6))
  })
  
  expect_message({
    fitted_model1 <- emuFit(Y = Y,
                           X = X[, 1, drop = FALSE],
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
                           test_kj = data.frame(k = 1, j = 1:6))
  })
  
  expect_equal(fitted_model$coef[, 2:9], fitted_model1$coef[, 2:9])
  
})

test_that("emuFit has 'score_test_hyperparams' object and throws warnings when convergence isn't hit", {
  # check that warning is returned when estimation under the alternative doesn't converge
  expect_warning({
    fitted_model <- emuFit(Y = Y,
                           X = cbind(X, rnorm(nrow(X))),
                           verbose = FALSE,
                           B_null_tol = 1e-2,
                           tolerance = 0.01,
                           tau = 2,
                           return_wald_p = FALSE,
                           compute_cis = FALSE,
                           run_score_tests = FALSE, 
                           maxit = 1)
  })
  expect_false(fitted_model$estimation_converged)
  
  # check that warning is returned when estimation under the null doesn't converge
  suppressWarnings({
    fitted_model <- emuFit(Y = Y,
                           X = cbind(X, rnorm(nrow(X))),
                           verbose = FALSE,
                           B_null_tol = 1e-2,
                           tolerance = 0.01,
                           tau = 2,
                           return_wald_p = FALSE,
                           compute_cis = FALSE,
                           run_score_tests = TRUE, 
                           test_kj = data.frame(k = 1, j = 1:2),
                           maxit = 1,
                           inner_maxit = 1)
  })
  
  # check that fitted model contains score_test_hyperparams object
  expect_true("score_test_hyperparams" %in% names(fitted_model))
  
  # check that fitted model contains data frame of unconverged test_kj
  expect_type(fitted_model$null_estimation_unconverged, "list")
})

test_that("test that B_null_list object can be used and throws appropriate warnings when used incorrectly", {
  expect_warning({
    fitted_model <- emuFit(Y = Y,
                           X = X,
                           B_null_list = list(b),
                           verbose = FALSE,
                           B_null_tol = 1e-2,
                           tolerance = 0.01,
                           tau = 2,
                           return_wald_p = FALSE,
                           compute_cis = FALSE,
                           run_score_tests = TRUE, 
                           test_kj = data.frame(k = 1, j = 1:2))
  })
  
  expect_silent({
    fitted_model <- emuFit(Y = Y,
                           X = X,
                           B_null_list = list(NULL, b),
                           verbose = FALSE,
                           B_null_tol = 1e-2,
                           tolerance = 0.01,
                           tau = 2,
                           return_wald_p = FALSE,
                           compute_cis = FALSE,
                           run_score_tests = TRUE, 
                           test_kj = data.frame(k = 1, j = 1:2))
  })
  
  expect_warning({
    fitted_model <- emuFit(Y = Y,
                           X = X,
                           B_null_list = list(NULL, b[, -3]),
                           verbose = FALSE,
                           B_null_tol = 1e-2,
                           tolerance = 0.01,
                           tau = 2,
                           return_wald_p = FALSE,
                           compute_cis = FALSE,
                           run_score_tests = TRUE, 
                           test_kj = data.frame(k = 1, j = 1:2))
  })
  
})

test_that("emuFit reorders X and X and Y rownames don't match", {
  dat1 <- data.frame(group = c(covariates$group[12], covariates$group[1:11]))
  rownames(dat1) <- paste0("sample", c(12, 1:11)) 
  dat2 <- covariates
  rownames(dat2) <- paste0("sample", 1:12)
  rownames(Y) <- paste0("sample", 1:12)
  
  expect_message({
    fitted_model1 <- emuFit(formula = ~ group,
                            data = dat1,
                            Y = Y,
                            compute_cis = FALSE,
                            run_score_tests = FALSE)
  })
  
  expect_silent({
    fitted_model2 <- emuFit(formula = ~ group,
                            data = dat2,
                            Y = Y,
                            compute_cis = FALSE,
                            run_score_tests = FALSE)
  })
  
  expect_true(all.equal(fitted_model1$coef, fitted_model2$coef))
})

test_that("emuFit throws error when there is a category with all zero counts", {
  
  Y_zero <- Y
  Y_zero[, 1] <- 0
  
  expect_error({
    fitted_model <- emuFit(Y = Y_zero,
                           X = X,
                           formula = ~group,
                           data = covariates,
                           verbose = FALSE,
                           B_null_tol = 1e-2,
                           tolerance = 0.01,
                           tau = 2,
                           return_wald_p = FALSE,
                           compute_cis = FALSE,
                           run_score_tests = FALSE)
  }) 
})

test_that("Confirm zi is different when penalty is applied or not", {
  
  fit_penT <- emuFit(Y = Y,
                     X = X,
                     formula = ~group,
                     data = covariates,
                     verbose = FALSE,
                     penalize = TRUE,
                     B_null_tol = 1e-2,
                     tolerance = 0.01,
                     tau = 2,
                     return_wald_p = FALSE,
                     compute_cis = FALSE,
                     run_score_tests = FALSE)
  
  fit_penF <- emuFit(Y = Y,
                     X = X,
                     formula = ~group,
                     data = covariates,
                     verbose = FALSE,
                     penalize = FALSE,
                     B_null_tol = 1e-2,
                     tolerance = 0.01,
                     tau = 2,
                     return_wald_p = FALSE,
                     compute_cis = FALSE,
                     run_score_tests = FALSE)
  
  expect_false(isTRUE(all.equal(fit_penT$z_hat,
                                fit_penF$z_hat)))
})

test_that("Single category constraint works", {
  
  emuRes <- emuFit(Y = Y,
                   X = X,
                   constraint_fn = 3,
                   run_score_tests = FALSE)
  expect_true(emuRes$B[2, 3] == 0)
  
})

test_that("emuFit works with fitted objects passed in", {
  emuRes <- emuFit(Y = Y,
                   X = X,
                   run_score_tests = FALSE)
  # can run emuFit with fitted model
  expect_silent({
    emuRes2 <- emuFit(Y = Y, X = X, fitted_model = emuRes, refit = FALSE,
                      compute_cis = FALSE, test_kj = data.frame(k = 2, j = 1))
  })
  # get error if have penalize arguments that don't match 
  expect_error({
    emuRes2 <- emuFit(Y = Y, X = X, fitted_model = emuRes, refit = FALSE,
                      compute_cis = FALSE, test_kj = data.frame(k = 2, j = 1),
                      penalize = FALSE)
  })
  # can run emuFit with only B 
  # can run emuFit with fitted model
  expect_silent({
    emuRes2 <- emuFit(Y = Y, X = X, B = emuRes$B, refit = FALSE,
                      compute_cis = FALSE, test_kj = data.frame(k = 2, j = 1))
  })
  
})

test_that("emuFit produces appropriate output when verbose = TRUE", {
  
  # check that when running score tests it tells you what is running and how long it has taken
  messages <- capture_messages(fitted_model <- emuFit(Y = Y,
                                                      X = X,
                                                      formula = ~group,
                                                      data = covariates,
                                                      verbose = TRUE,
                                                      B_null_tol = 1e-2,
                                                      tolerance = 0.01,
                                                      tau = 2,
                                                      test_kj = data.frame(j = 1, k = 2)))
  expect_true(TRUE %in% grepl("Running score test", messages) & 
                TRUE %in% grepl("has completed in approximately", messages))
  
})

test_that("emuFit refits starting at provided value if `B` or `fitted_model` are given", {
  
  message1 <- capture_messages(fitted_model <- emuFit(Y = Y,
                                                      X = X,
                                                      formula = ~group,
                                                      data = covariates,
                                                      B_null_tol = 1e-2,
                                                      tolerance = 0.01,
                                                      tau = 2,
                                                      run_score_tests = FALSE,
                                                      compute_cis = FALSE, 
                                                      verbose = "development"))
  message2 <- capture_messages(fitted_model_with_B <- emuFit(Y = Y,
                                                             X = X,
                                                             B = fitted_model$B, 
                                                             refit = TRUE, 
                                                             formula = ~group,
                                                             data = covariates,
                                                             B_null_tol = 1e-2,
                                                             tolerance = 0.01,
                                                             tau = 2,
                                                             run_score_tests = FALSE,
                                                             compute_cis = FALSE,
                                                             verbose = "development"))
  message3 <- capture_messages(fitted_model_with_fit <- emuFit(Y = Y,
                                                               X = X,
                                                               fitted_model = fitted_model, 
                                                               refit = TRUE,
                                                               formula = ~group,
                                                               data = covariates,
                                                               B_null_tol = 1e-2,
                                                               tolerance = 0.01,
                                                               tau = 2,
                                                               run_score_tests = FALSE,
                                                               compute_cis = FALSE,
                                                               verbose = "development"))
  expect_true(length(message1) > length(message2))
  expect_true(length(message1) > length(message3))
  
})

test_that("giving test_kj as valid strings works", {
  colnames(Y) <- paste0("taxon", 1:6)
  colnames(X) <- c("int", "group")
  res <- emuFit(Y = Y, X = X, compute_cis = FALSE, test_kj = data.frame(k = "group", j = "taxon3"),
                penalize = FALSE, tolerance = 0.1)
  expect_true(!is.na(res$coef$pval[3]))
  expect_error(emuFit(Y = Y, X = X, compute_cis = FALSE, test_kj = data.frame(k = "group", j = "taxa3"),
                      penalize = FALSE, tolerance = 0.1))
})
