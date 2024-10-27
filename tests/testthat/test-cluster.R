test_that("clusters work as I want", {
  
  set.seed(100)
  n <- 64
  J <- 50
  
  cage_num <- rep(c(1:16), 4) 
  treatment <- (cage_num <= 8)
  XX <- data.frame(treatment)
  
  Y <- radEmu:::simulate_data(n = n,
                              J = J,
                              X = cbind(1, treatment),
                              b0 = runif(J, min = 0, max = 4),
                              b1 = runif(J, min = 0, max = 4),
                              distn = "ZINB",
                              zinb_size = 10,
                              zinb_zero_prop = 0.3,
                              mean_z = 5)
  
  # check that cluster argument works as a numeric vector 
  ef_num <- emuFit(formula = ~ treatment, 
                   data = XX, 
                   Y = Y, 
                   cluster=cage_num, 
                   run_score_tests=FALSE) #### very fast
  expect_equal(ef_num$coef %>% class, "data.frame")
  
  # check that cluster argument works as character vector and gives 
  # equivalent results to numeric vector 
  cage_char <- rep(c(LETTERS[1:16]), 4)
  ef_char <- emuFit(formula = ~ treatment, 
                    data = XX, 
                    Y = Y, 
                    cluster=cage_char, 
                    run_score_tests=FALSE) 
  expect_equal(ef_num$coef, ef_char$coef)
  
  # check that cluster argument works as factor and gives equivalent results
  # to numeric vector 
  cage_fact <- cage_char %>% as.factor
  ef_fact <- emuFit(formula = ~ treatment, 
                    data = XX, 
                    Y = Y, 
                    cluster=cage_fact, 
                    run_score_tests=FALSE) 
  expect_equal(ef_num$coef, ef_fact$coef)
})


set.seed(11)
X <- cbind(1,rnorm(12))
Y <- radEmu:::simulate_data(n = 12,
                            J = 6,
                            X = X,
                            b0 = rnorm(J),
                            b1 = seq(1,5,length.out = 6) - mean(seq(1,5,length.out = 6)),
                            distn = "ZINB",
                            zinb_size = 2,
                            zinb_zero_prop = 0.8,
                            mean_z = 5)
covariates <- data.frame(group = X[,2])


test_that("GEE with cluster covariance gives plausible type 1 error ",{
  
  skip("Skipping -- a simulation for T1E under cluster dependence.")
  
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
    # print(sim)
    
    X <- cbind(1, rnorm(12))
    Y <- radEmu:::simulate_data(n = 12,
                                J = 6,
                                X = X,
                                b0 = rnorm(J),
                                b1 = seq(1,5,length.out = J) - mean(seq(1,5,length.out = J)),
                                distn = "ZINB",
                                zinb_size = 5,
                                zinb_zero_prop = 0.8,
                                mean_z = 5,
                                cluster = cluster)
    covariates <- data.frame(group = X[,2])
    
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
  
  #expect conservative inference for score test with correct clustering
  expect_true(mean(results$pval<=0.05) <= 0.05)
  #expect anti-conservative inference for score test without clustering
  expect_true(mean(results_noGEE$pval<=0.05) > 0.06)
  
})
