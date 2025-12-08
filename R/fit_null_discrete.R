fit_null_discrete_pseudohuber <- function(
    B,
    Y,
    X,
    k_constr,
    j_constr,
    j_ref,
    maxit = 1000
) {
  
  # if (is.null(j_ref)) {
  #   j_ref <- ifelse(j_constr == 1, 2, 1)
  # }
  
  n <- nrow(Y)
  J <- ncol(Y)
  p <- ncol(X)
  
  ## TODO generalise
  stopifnot(j_ref == J)
  
  distinct_xx <- unique(X)
  
  stopifnot(ncol(X) == nrow(distinct_xx))
  stopifnot(all(X[,1] == 1)) ### need a baseline category to simplify the parametrization
  X_wo_1s <- distinct_xx[, 1]
  
  groups <- split(
    seq_len(nrow(X)),                        # row indices of original data
    apply(X, 1, function(r)
      paste(r, collapse = "_"))               # grouping key by row contents
  )
  
  groups <- groups[order(sapply(groups, min))]
  
  totals <- lapply(groups, function(x) { 
    apply(Y[x, , drop = FALSE], 2, sum)
  })
  
  Y_sum <- do.call(rbind, totals)
  
  ## this should be equivalent to fit_null with j_ref=J
  out <- fit_null_discrete_micro_fs(
    Y = Y_sum, 
    X = X_wo_1s, 
    k_constr=k_constr, 
    j_constr=j_constr,
    constraint_fn=function(x) {  pseudohuber_median(c(x, 0)) },       
    constraint_grad_fn= function(x) {  x <- radEmu::dpseudohuber_median_dx(c(x, 0)); x[-length(x)]}
    maxit = maxit,
  )
  
  ## Currently this works for 
  #### two groups with X's (1, 0) and (1, 1)
  #### testing the first column, last column also constrained
  
  #### it is equivalent to fit_null with k_constr=2, j_constr=1, j_ref=J,
  
  ## It needs to be generalized to consider 
  #### multiple categories (can we assume identity-like then back transform? is this constraint and likelihood-preserving?) 
  #### different columns for testing
  #### the same inputs as fit_null and comparable convergence statistics
  
  # out <- my_fs_stable_two_groups(n0=Y[which(X[,2] == 0), ] %>% colSums, 
  #                                n1 = Y[which(X[,2] == 1), ] %>% colSums, 
  #                                g_beta=function(x) {  pseudohuber_median(c(x, 0)) },  
  #                                g_beta_grad= function(x) {  x <- radEmu::dpseudohuber_median_dx(c(x, 0)); x[-length(x)]}, 
  #                                maxit = 1000, tol = 1e-8,
  #                                ls_max = 20, ls_rho = 0.5,
  #                                ridge_base = 1e-4,
  #                                max_step_norm = 5,
  #                                clip_logit_max = 15)
  
  B <- out$B
  
  z <- update_z(Y, X, B)
  log_means <- X %*% B + matrix(z, ncol = 1) %*% matrix(1, ncol = J, nrow = 1)
  ll_new <- sum(Y * log_means - exp(log_means))
  
  # beta = X^{-1} %*% eta
  betahats <- round(MASS::ginv(distinct_xx), 8) %*% etahats
  betahats
  
  
  it_df <- it_df[1:iter, ]
  
  return(list(
    "B" = B,
    "k_constr" = k_constr,
    "j_constr" = j_constr,
    "niter" = iter,
    "gap" = gap,
    "u" = u,
    "rho" = rho,
    "Bs" = Bs,
    "it_df" = it_df,
    "converged" = converged
  ))
  
}
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
