test_that("clusters work as I want", {
  
  set.seed(100)
  n <- 64
  J <- 50
  Y <- matrix(rpois(n*J, 100)*rbinom(n*J, 100, 0.7), nrow=n, ncol = J)
  cage_num <- rep(c(1:16), 4) 
  treatment <- (cage_num <= 8)
  XX <- data.frame(treatment)
  
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
covariates <- data.frame(group = X[,2])

b0 <- rnorm(J)
b1 <- seq(1,5,length.out = J)
b1 <- b1 - mean(b1)
b1[3:4] <- 0
b <- rbind(b0,b1)


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
  
  #expect conservative inference for score test with correct clustering
  expect_true(mean(results$pval<=0.05) <= 0.05)
  #expect anti-conservative inference for score test without clustering
  expect_true(mean(results_noGEE$pval<=0.05) > 0.06)
  
})
