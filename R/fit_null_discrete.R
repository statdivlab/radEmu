# fit_null_discrete <- function(
#     B,
#     Y,
#     X,
#     k_constr,
#     j_constr,
#     j_ref,
#     constraint_fn,
#     constraint_grad_fn
# ) {
#   
#   if (is.null(j_ref)) {
#     j_ref <- ifelse(j_constr == 1, 2, 1)
#   }
#   
#   n <- nrow(Y)
#   J <- ncol(Y)
#   p <- ncol(X)
#   
#   distinct_xx <- unique(X)
#   
#   # fitted values eta = X %*% beta
#   pihats <- matrix(NA, nrow = p, ncol = J) # p x J
#   etahats <- matrix(NA, nrow = p, ncol = J) # p x J
#   etahats[, j_ref] <- 0
#   
#   groups <- split(
#     seq_len(nrow(X)),                        # row indices of original data
#     apply(X, 1, function(r)
#       paste(r, collapse = "_"))               # grouping key by row contents
#   )
#   groups <- groups[order(sapply(groups, min))]
#   
#   for (the_cat in 1:length(groups)) {
#     the_xs <- groups[[the_cat]]
#     totals <- apply(Y[the_xs, , drop = FALSE], 2, sum)
#     pihats[the_cat, ] <- totals / sum(totals)
#     etahats[the_cat, setdiff(1:J, j_ref)] <- log(pihats[-j_ref] / pihats[j_ref])
#   }
#   
#   
#   
#   
#   # beta = X^{-1} %*% eta
#   betahats <- round(MASS::ginv(distinct_xx), 8) %*% etahats
#   betahats
#   
# }
# 
# 
# 
# set.seed(1)
# n <- 10
# J <- 5
# p <- 2
# beta <- matrix(rnorm(p*J, sd=5), ncol = J)
# X <- cbind(1, c(rep(0, n/2), rep(1, n/2)))
# Y <- matrix(rpois(n=n*J, lambda=exp(5 + X %*% beta)), ncol = J)
# Y
# 
# # woah, that takes ages
# fit1 <- radEmu::emuFit(Y = Y, X = X, test_kj=data.frame("k" = 2, "j" = 3))
# fit1
# fit1 %>% names
# fit1$estimation_converged
# fit1$score_test_hyperparams
# fit1
# 
# fit2 <- radEmu::emuFit(Y = Y, X = X, test_kj=data.frame("k" = 2, "j" = 3), control=list("trackB" = TRUE), refit=F, fitted_model=fit1)
# fit2
# fit2 %>% names
# library(tidyverse)
# final_bs <- fit2$trackB_list[[1]] %>% tibble %>% tail(p*J)
# final_bs %>% 
#   signif(3)
# fit2$coef 
# 
# 
# final_bs %>% 
#   signif(3)
# Y
# X
# final_bs %>% pull(B) %>% matrix(nrow = 2, byrow=T)
# 
# 
# #### attempt 2
# set.seed(1)
# n <- 10
# J <- 5
# p <- 2
# beta <- matrix(rnorm(p*J, sd=2), ncol = J)
# X <- cbind(1, c(rep(0, n/2), rep(1, n/2)))
# Y <- matrix(rpois(n=n*J, lambda=exp(5 + X %*% beta)), ncol = J)
# Y
# 
# fit3 <- radEmu::emuFit(Y = Y, X = X, test_kj=data.frame("k" = 2, "j" = 4), control=list("trackB" = TRUE))
# fit3
# fit3$trackB_list[[1]] %>% tibble %>% tail(p*J)
# 
# fit3$trackB_list[[1]] %>% tibble %>% tail(J) %>% pull(B) %>% psuedohuber_median
# fit3$trackB_list[[1]] %>% tibble %>% tail(p*J) %>% head(J) %>% pull(B) %>% psuedohuber_median
# 
# final_bs %>% 
#   signif(3)
# Y
# X
# final_bs %>% pull(B) %>% matrix(nrow = 2, byrow=T)
# 
# 
# fit3$B
# 
# fn3b <- fit_null(B=fit3$B, Y=Y, X = X, k_constr=2, j_constr=4, j_ref=2, 
#                  constraint_fn=list(pseudohuber_median, pseudohuber_median), 
#                  constraint_grad_fn=list(radEmu::dpseudohuber_median_dx, radEmu::dpseudohuber_median_dx))
# fn3a$B
# fn3b$B
# 
# fn3a <- fit_null(B=fit3$B, Y=Y, X = X, k_constr=2, j_constr=4, j_ref=1, 
#                  constraint_fn=list(pseudohuber_median, pseudohuber_median), 
#                  constraint_grad_fn=list(radEmu::dpseudohuber_median_dx, radEmu::dpseudohuber_median_dx), 
#                  B_tol=1e-8, constraint_tol=1e-8)
# totest <- fn3a
# totest$B
# totest$B[totest$k_constr, setdiff(1:J, totest$j_constr)] %>% psuedohuber_median
# totest$B[totest$k_constr, totest$j_constr] 
# 
# 
# 
# aaa <- fit3$trackB_list[[1]] %>% tibble %>% tail(p*J) %>% tail(J)
# (aaa[c(1,2,3,5)] %>% psuedohuber_median) - aaa[4] 
