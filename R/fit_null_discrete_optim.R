## n0: counts for x = 0, length J
## n1: counts for x = 1, length J
## gfun: function(beta_free2) -> scalar beta_1
##       beta_free2 is c(beta_2, ..., beta_{J-1})
## gprimefun: function(beta_free2) -> numeric(J-2)
##            gradient d beta_1 / d (beta_2,...,beta_{J-1})

loglik_counts_constrained <- function(par, n0, n1, gfun, gprimefun) {
  J  <- length(n0)
  Jm1 <- J - 1
  
  ## Split parameters:
  ## alpha_1,...,alpha_{J-1}, beta_2,...,beta_{J-1}
  alpha       <- par[1:Jm1]
  beta_free2  <- par[(J):(2*J - 3)]   # length J-2: beta_2..beta_{J-1}
  
  ## Constrained beta_1 via g
  beta1 <- gfun(beta_free2)
  
  ## Full beta for categories 1..J-1 (baseline J has beta_J = 0)
  beta_full <- numeric(Jm1)
  beta_full[1] <- beta1
  beta_full[2:Jm1] <- beta_free2
  
  ## Group sizes
  N0 <- sum(n0)
  N1 <- sum(n1)
  
  ## x = 0: eta0_j = alpha_j, eta0_J = 0
  eta0 <- alpha                       # length J-1
  m0   <- max(c(0, eta0))
  exp_eta0 <- exp(eta0 - m0)
  denom0   <- exp(0 - m0) + sum(exp_eta0)
  p0       <- exp_eta0 / denom0      # probs for j=1..J-1
  p0J      <- exp(0 - m0) / denom0   # prob for baseline J
  
  ## x = 1: eta1_j = alpha_j + beta_j, eta1_J = 0
  eta1 <- alpha + beta_full
  m1   <- max(c(0, eta1))
  exp_eta1 <- exp(eta1 - m1)
  denom1   <- exp(0 - m1) + sum(exp_eta1)
  p1       <- exp_eta1 / denom1      # probs for j=1..J-1
  p1J      <- exp(0 - m1) / denom1
  
  ## Log-likelihood:
  ## ℓ = sum_j n0_j log p0_j + sum_j n1_j log p1_j
  ## with p0_J, p1_J included for category J
  ll0 <- sum(n0[1:Jm1] * log(p0))  + n0[J] * log(p0J)
  ll1 <- sum(n1[1:Jm1] * log(p1))  + n1[J] * log(p1J)
  
  ll0 + ll1
}

grad_counts_constrained <- function(par, n0, n1, gfun, gprimefun) {
  J  <- length(n0)
  Jm1 <- J - 1
  
  alpha      <- par[1:Jm1]
  beta_free2 <- par[(J):(2*J - 3)]
  
  beta1 <- gfun(beta_free2)
  beta_full <- numeric(Jm1)
  beta_full[1] <- beta1
  beta_full[2:Jm1] <- beta_free2
  
  gprime <- gprimefun(beta_free2)   # length J-2
  
  N0 <- sum(n0)
  N1 <- sum(n1)
  
  ## x = 0
  eta0 <- alpha
  m0   <- max(c(0, eta0))
  exp_eta0 <- exp(eta0 - m0)
  denom0 <- exp(0 - m0) + sum(exp_eta0)
  p0 <- exp_eta0 / denom0                  # j = 1..J-1
  p0J <- exp(0 - m0) / denom0
  
  ## x = 1
  eta1 <- alpha + beta_full
  m1   <- max(c(0, eta1))
  exp_eta1 <- exp(eta1 - m1)
  denom1 <- exp(0 - m1) + sum(exp_eta1)
  p1 <- exp_eta1 / denom1                  # j = 1..J-1
  p1J <- exp(0 - m1) / denom1
  
  ## Gradient w.r.t. alpha_j (j=1..J-1):
  ## ∂ℓ/∂alpha_j = (n0_j - N0 p0_j) + (n1_j - N1 p1_j)
  grad_alpha <- (n0[1:Jm1] - N0 * p0) + (n1[1:Jm1] - N1 * p1)
  
  ## Unconstrained slope scores:
  ## ∂ℓ/∂beta_j = n1_j - N1 p1_j, j = 1..J-1
  score_beta_full <- n1[1:Jm1] - N1 * p1
  
  ## ∂ℓ/∂beta_1:
  S1 <- score_beta_full[1]
  
  ## For free beta_2..beta_{J-1}:
  ## ∂ℓ/∂beta_j_free = (n1_j - N1 p1_j) + S1 * ∂beta1/∂beta_j
  ## j = 2..J-1
  S_direct <- score_beta_full[2:Jm1]       # length J-2
  grad_beta_free2 <- S_direct + S1 * gprime
  
  c(grad_alpha, grad_beta_free2)
}

fit_multinom_constrained_counts <- function(n0, n1, gfun, gprimefun,
                                            start = NULL) {
  J  <- length(n0)
  if (length(n1) != J) stop("n0 and n1 must have the same length")
  if (J < 3) stop("Need at least 3 categories (J >= 3).")
  
  Jm1 <- J - 1
  n_par <- (Jm1) + (Jm1 - 1)   # alpha_1..alpha_{J-1}, beta_2..beta_{J-1}
  
  if (is.null(start)) {
    start <- rep(0, n_par)
  }
  
  fn <- function(par) {
    -loglik_counts_constrained(par, n0, n1, gfun, gprimefun)  # optim minimizes
  }
  
  gr <- function(par) {
    -grad_counts_constrained(par, n0, n1, gfun, gprimefun)
  }
  
  opt <- optim(start, fn, gr, method = "BFGS", control = list(maxit = 100))
  
  par_hat <- opt$par
  
  ## Reconstruct alpha_1..alpha_J and beta_1..beta_J
  alpha_hat <- par_hat[1:Jm1]
  beta_free2_hat <- par_hat[(J):(2*J - 3)]
  beta1_hat <- gfun(beta_free2_hat)
  
  beta_hat <- numeric(J)
  beta_hat[1] <- beta1_hat
  beta_hat[2:Jm1] <- beta_free2_hat
  beta_hat[J] <- 0
  
  alpha_hat_full <- numeric(J)
  alpha_hat_full[1:Jm1] <- alpha_hat
  alpha_hat_full[J] <- 0
  
  list(
    alpha = alpha_hat_full,
    beta  = beta_hat,
    logLik = -opt$value,
    convergence = opt$convergence,
    optim = opt
  )
}

# system.time({
# out_test3 <- 
fit_multinom_constrained_counts(n0=Y[which(X[,2] == 0), ] %>% colSums, 
                                n1 = Y[which(X[,2] == 1), ] %>% colSums, 
                                gfun=function(x) {  pseudohuber_median(c(x, 0)) },
                                gprimefun= function(x) {  x <- radEmu::dpseudohuber_median_dx(c(x, 0)); x[-length(x)]}
)
# })
out_test2$alpha
out_test2$beta
out_test3
