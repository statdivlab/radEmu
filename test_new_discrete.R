# js <- seq(from = 100, to = 1000, by = 100)
js <- seq(from = 20, to = 100, by = 20)
df_timing <- data.frame("j" = js, "fit_null" = NA, "new" = NA)
df_timing
set.seed(1)
ii <- 1
for (jj in js) {
  df_timing$j[ii] <- jj
  
  n <- 20
  J <- jj
  p <- 2
  beta <- matrix(rnorm(p*J, sd=5), ncol = J)
  X <- cbind(1, c(rep(0, n/2), rep(1, n/2)))
  Y <- matrix(rpois(n=n*J, lambda=exp(5 + X %*% beta)), ncol = J)
  Y <- pmax(Y, 1)
  
  # fit3 <- radEmu::emuFit(Y = Y, X = X, run_score_tests=F)
  # 
  # df_timing$fit_null[ii] <- system.time({
  #   fit_null(B=fit3$B, Y=Y, X = X, k_constr=2, j_constr=1, j_ref=J,
  #            constraint_fn=list(pseudohuber_median, pseudohuber_median),
  #            constraint_grad_fn=list(radEmu::dpseudohuber_median_dx, radEmu::dpseudohuber_median_dx),
  #            B_tol=1e-8, constraint_tol=1e-8)
  # })[3]
  
  
  df_timing$fit_null[ii] <- system.time({
    fn <- fit_null(B=matrix(0, nrow = 2, ncol = J), Y=Y, X = X, k_constr=2, j_constr=1, j_ref=J,
                   constraint_fn=list(pseudohuber_median, pseudohuber_median),
                   constraint_grad_fn=list(radEmu::dpseudohuber_median_dx, radEmu::dpseudohuber_median_dx),
                   B_tol=1e-5, constraint_tol=1e-8)
  })[3]
  
  # df_timing$new[ii] <- system.time({
  #   out_test4 <- my_gd_fs(n0=Y[which(X[,2] == 0), ] %>% colSums, 
  #                         n1 = Y[which(X[,2] == 1), ] %>% colSums, 
  #                         g_beta=function(x) {  pseudohuber_median(c(x, 0)) },  
  #                         g_beta_grad= function(x) {  x <- radEmu::dpseudohuber_median_dx(c(x, 0)); x[-length(x)]}, 
  #                         eta_alpha = 1e-3,
  #                         eta_beta  = 1e-3,
  #                         maxit = 1e10,
  #                         B_tol = 1e-5)
  # })[3]
  df_timing$new[ii] <- system.time({
    out_test4 <- my_fs_stable_two_groups(n0=Y[which(X[,2] == 0), ] %>% colSums, 
                                         n1 = Y[which(X[,2] == 1), ] %>% colSums, 
                                         g_beta=function(x) {  pseudohuber_median(c(x, 0)) },  
                                         g_beta_grad= function(x) {  x <- radEmu::dpseudohuber_median_dx(c(x, 0)); x[-length(x)]}, 
                                         eta_alpha = 1e-3,
                                         eta_beta  = 1e-3,
                                         maxit = 1e8,
                                         tol = 1e-10)
  })[3]
  ii <- ii + 1
  print(df_timing)
}
df_timing

out_test4
rbind(out_test4$alpha, out_test4$beta)
fn$B
### why didn't they align???


set.seed(1)
n <- 20
J <- 100
p <- 2
beta <- matrix(rnorm(p*J, mean = 2, sd=2), ncol = J)
X <- cbind(1, c(rep(0, n/2), rep(1, n/2)))
Y <- matrix(rpois(n=n*J, lambda=exp(2 + X %*% beta)), ncol = J)
Y <- pmax(Y, 1)
Y
Ysimple <- rbind(Y[which(X[,2] == 0), ] %>% colSums, Y[which(X[,2] == 1), ] %>% colSums)

# the_fit <- radEmu::emuFit(Y = Y, X = X, run_score_tests=F)

the_fit <- radEmu::emuFit(Y = Ysimple, X = X[10:11, ], run_score_tests=F)
the_fit

system.time({
  fn3_1_J <- fit_null(B=the_fit$B, Y=Y, X = X, k_constr=2, j_constr=1, j_ref=J,
                      constraint_fn=list(pseudohuber_median, pseudohuber_median),
                      constraint_grad_fn=list(radEmu::dpseudohuber_median_dx, radEmu::dpseudohuber_median_dx),
                      B_tol=1e-8, constraint_tol=1e-8)
}) #  26.270   0.152  27.411 

system.time({
  out_test2 <- my_gd_ls(n0=Y[which(X[,2] == 0), ] %>% colSums,
                        n1 = Y[which(X[,2] == 1), ] %>% colSums,
                        g_beta=function(x) {  pseudohuber_median(c(x, 0)) },
                        g_beta_grad= function(x) {  x <- radEmu::dpseudohuber_median_dx(c(x, 0)); x[-length(x)]},
                        eta_alpha = 1e-3,
                        eta_beta  = 1e-3,
                        maxit = 1000,
                        tol = 1e-6)
}) # 2.912   0.033   3.287 


system.time({
  out_test4 <- my_gd_fs(n0=Y[which(X[,2] == 0), ] %>% colSums,
                        n1 = Y[which(X[,2] == 1), ] %>% colSums,
                        g_beta=function(x) {  pseudohuber_median(c(x, 0)) },
                        g_beta_grad= function(x) {  x <- radEmu::dpseudohuber_median_dx(c(x, 0)); x[-length(x)]},
                        eta_alpha = 1e-3,
                        eta_beta  = 1e-3,
                        maxit = 1000,
                        B_tol = 1e-8)
}) #  0.021   0.001   0.023 


system.time({
  out_test5 <- my_fs_stable(n0=Y[which(X[,2] == 0), ] %>% colSums,
                            n1 = Y[which(X[,2] == 1), ] %>% colSums,
                            g_beta=function(x) {  pseudohuber_median(c(x, 0)) },
                            g_beta_grad= function(x) {  x <- radEmu::dpseudohuber_median_dx(c(x, 0)); x[-length(x)]},
                            eta_alpha = 1e-3,
                            eta_beta  = 1e-3,
                            maxit = 1000,
                            tol = 1e-8)
}) # 0.016   0.000   0.017 



# out_test4$history %>% filter(type == "beta" & k == 3)
# 
fn3_1_J$B
rbind(out_test2$alpha, out_test2$beta)
rbind(out_test4$alpha, out_test4$beta)
rbind(out_test5$alpha, out_test5$beta)

####################################################
#### Try larger sample


set.seed(1)
n <- 20
J <- 100
p <- 2
beta <- matrix(rnorm(p*J, mean = 2, sd=2), ncol = J)
X <- cbind(1, c(rep(0, n/2), rep(1, n/2)))
Y <- matrix(rpois(n=n*J, lambda=exp(2 + X %*% beta)), ncol = J)
Y <- pmax(Y, 1)
Y
Ysimple <- rbind(Y[which(X[,2] == 0), ] %>% colSums, Y[which(X[,2] == 1), ] %>% colSums)

system.time({
  fn3_1_J <- fit_null(B=matrix(0, nrow = p, ncol = J), Y=Y, X = X, k_constr=2, j_constr=1, j_ref=J,
                      constraint_fn=list(pseudohuber_median, pseudohuber_median),
                      constraint_grad_fn=list(radEmu::dpseudohuber_median_dx, radEmu::dpseudohuber_median_dx),
                      B_tol=1e-8, constraint_tol=1e-8)
}) 

system.time({
  out_test5 <- my_fs_stable_two_groups(n0=Y[which(X[,2] == 0), ] %>% colSums,
                                       n1 = Y[which(X[,2] == 1), ] %>% colSums,
                                       g_beta=function(x) {  pseudohuber_median(c(x, 0)) },
                                       g_beta_grad= function(x) {  x <- radEmu::dpseudohuber_median_dx(c(x, 0)); x[-length(x)]},
                                       maxit = 1000,
                                       tol = 1e-8)
}) 

orig <- fn3_1_J$B
new <- cbind(rbind(out_test5$alpha, out_test5$beta), 0)

## are the constraints satisfied
c(psuedohuber_median(orig[2,-1]), psuedohuber_median(orig[2,1]))
c(psuedohuber_median(new[2,-1]), psuedohuber_median(new[2,1]))

B <- orig
z <- update_z(Y, X, B)
log_means <- X %*% B + matrix(z, ncol = 1) %*% matrix(1, ncol = J, nrow = 1)
ll_orig <- sum(Y * log_means - exp(log_means))

B <- new
z <- update_z(Y, X, B)
log_means <- X %*% B + matrix(z, ncol = 1) %*% matrix(1, ncol = J, nrow = 1)
ll_new <- sum(Y * log_means - exp(log_means))

B <- matrix(0, nrow = p, ncol = J)
z <- update_z(Y, X, B)
log_means <- X %*% B + matrix(z, ncol = 1) %*% matrix(1, ncol = J, nrow = 1)
sum(Y * log_means - exp(log_means))

ll_orig < ll_new ### new is better


##### try g-dim version

set.seed(1)
n <- 20
J <- 30
p <- 3
beta <- matrix(rnorm(p*J, mean = 2, sd=1), ncol = J)
X <- cbind(1, c(rep(0, n/2), rep(1, n/2)), c(rep(0, 3*n/4), rep(1, n/4)))
Y <- matrix(rpois(n=n*J, lambda=exp(2 + X %*% beta)), ncol = J)
Y <- pmax(Y, 1)
Y

# radEmu:::emuFit_micro_discrete(X, Y)
# radEmu:::emuFit_micro_discrete(X[n:1, ], Y[n:1, ])

unique(X)

my_n_list <- list(Y[apply(X, 1, function(x) all(x == unique(X)[1,])), ] %>% colSums, 
                  Y[apply(X, 1, function(x) all(x == unique(X)[2,])), ] %>% colSums, 
                  Y[apply(X, 1, function(x) all(x == unique(X)[3,])), ] %>% colSums)


load_all()
my_jstar = 9
system.time({
  out_test6 <- fs_direct_beta_constraint(
    Y = do.call(rbind, my_n_list), 
    X = unique(X)[, -1], k_constr=2, j_constr=my_jstar,
    constraint_fn=function(x) {  pseudohuber_median(c(x, 0)) },       # f(z2..zm)
    constraint_grad_fn= function(x) {  x <- radEmu::dpseudohuber_median_dx(c(x, 0)); x[-length(x)]}
  )
}) 
out_test6

out_test6$iter
out_test6$alpha %>% c(0) %>% pseudohuber_median

(out_test6$beta[,-my_jstar] %>% apply(1, function(x) {  pseudohuber_median(c(x, 0)) }))[2]
out_test6$beta[2,my_jstar] 
### ok, promising



system.time({
  out_test6_orig <- fit_null(B=matrix(0, nrow = p, ncol = J), 
                             Y=Y, X = X, 
                             k_constr=2, j_constr=my_jstar, j_ref=J,
                             constraint_fn=list(pseudohuber_median, pseudohuber_median),
                             constraint_grad_fn=list(radEmu::dpseudohuber_median_dx, radEmu::dpseudohuber_median_dx),
                             B_tol=1e-8, constraint_tol=1e-8)
}) # 246.458   0.906 247.997
out_test6_orig$B
cbind(rbind(out_test6$alpha, out_test6$beta), 0)
## excellent 
## now check constraints
