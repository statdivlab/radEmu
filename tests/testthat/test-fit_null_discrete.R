
test_that("new discrete is correct", {
  
  set.seed(1)
  n <- 20
  J <- 10
  p <- 3
  beta <- matrix(rnorm(p*J, mean = 2, sd=1), ncol = J)
  X <- cbind(1, c(rep(0, n/2), rep(1, n/2)), c(rep(0, 3*n/4), rep(1, n/4)))
  Y <- matrix(rpois(n=n*J, lambda=exp(2 + X %*% beta)), ncol = J)
  Y <- pmax(Y, 1)
  
  my_n_list <- list(Y[apply(X, 1, function(x) all(x == unique(X)[1,])), ] %>% colSums, 
                    Y[apply(X, 1, function(x) all(x == unique(X)[2,])), ] %>% colSums, 
                    Y[apply(X, 1, function(x) all(x == unique(X)[3,])), ] %>% colSums)
  
  
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
  
  
  
})