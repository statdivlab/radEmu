# test_that("robust score p-values are reasonable in marginal fit", {
#
#   n <- 50
#   X <- cbind(1,c(rep(0,n),rep(1,n)))
#   set.seed(43234)
#   z <- rnorm(2*n,7,3)
#   log_means <- cbind(z,z + 3)
#   pvals <- numeric(100)
#   score_stats <- numeric(100)
#   for(i in 1:100){
#     print(i)
#     Y <- apply(log_means,c(1,2),
#                function(x) rpois(1,exp(x)))
#                  # rnbinom(1,size = 2, mu = exp(x)))
#     obj <- try(fit_marginal(X,Y))
#     pvals[i] <- try(obj$params$pval[2])
#     score_stats[i] <- try(obj$params$score_stat[2])
#   }
# })
