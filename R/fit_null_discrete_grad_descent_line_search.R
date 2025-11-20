my_gd_ls  <- function(n0, n1,
                      g_beta,  g_beta_grad,
                      eta_alpha = 1e-3,
                      eta_beta  = 1e-3,
                      maxit = 1000,
                      tol = 1e-8,
                      ls_max = 30,        # max backtracking steps
                      ls_c   = 1e-4,      # Armijo constant
                      ls_rho = 0.5) {     # backtracking shrink factor
  
  K  <- length(n0)
  Km <- K - 1
  
  if (K < 3) {
    stop("beta constraint beta1 = g(beta2,...,beta_{K-1}) requires K >= 3.")
  }
  
  # Parameters for categories 1..K-1 (baseline K is 0)
  alpha <- rep(0, Km)
  beta  <- rep(0, Km)
  
  # Enforce beta constraint initially: beta1 = g(beta2..beta_{K-1})
  beta[1] <- g_beta(beta[2:Km])
  
  n0p <- sum(n0)
  n1p <- sum(n1)
  
  # Softmax probabilities for x = 0 and x = 1 (excluding baseline)
  p0fun <- function(a) {
    ea <- exp(a)
    ea / (1 + sum(ea))
  }
  p1fun <- function(a, b) {
    eb <- exp(a + b)
    eb / (1 + sum(eb))
  }
  
  # Log-likelihood given alpha, beta
  loglik_fun <- function(alpha, beta) {
    p0 <- p0fun(alpha)
    p1 <- p1fun(alpha, beta)
    
    p0_full <- c(p0, 1 - sum(p0))
    p1_full <- c(p1, 1 - sum(p1))
    
    # clip probs to avoid log(0)
    eps <- 1e-12
    p0_full <- pmin(pmax(p0_full, eps), 1 - eps)
    p1_full <- pmin(pmax(p1_full, eps), 1 - eps)
    
    sum(n0 * log(p0_full)) + sum(n1 * log(p1_full))
  }
  
  # ---- history tracking ----
  param_history <- list()
  # iteration 0 (initial)
  param_history[[length(param_history) + 1]] <- rbind(
    data.frame(iter = 0L,
               type = "alpha",
               k = seq_len(Km),
               value = alpha),
    data.frame(iter = 0L,
               type = "beta",
               k = seq_len(Km),
               value = beta)
  )
  # ---------------------------
  
  # Initial log-likelihood
  ll_old <- loglik_fun(alpha, beta)
  
  final_iter <- 0L
  
  for (iter in 1:maxit) {
    
    # Probabilities for non-baseline categories 1..K-1
    p0 <- p0fun(alpha)
    p1 <- p1fun(alpha, beta)
    
    # Unconstrained gradients wrt alpha[1:Km] and beta[1:Km]
    dlda <- (n0[1:Km] + n1[1:Km]) - n0p * p0 - n1p * p1
    dldb <-  n1[1:Km]            - n1p * p1
    
    # Gradient for beta free coordinates k = 2..Km (chain rule)
    grad_g_beta     <- g_beta_grad(beta[2:Km])             # length K-2
    grad_beta_free  <- dldb[2:Km] + dldb[1] * grad_g_beta  # length K-2
    
    # Gradient norm for convergence check (use α and free β)
    grad_norm <- sqrt(sum(dlda^2) + sum(grad_beta_free^2))
    if (grad_norm < tol) {
      final_iter <- iter - 1L
      break
    }
    
    # Search directions (gradient ascent)
    dir_alpha     <- eta_alpha * dlda
    dir_beta_free <- eta_beta  * grad_beta_free
    
    # Directional magnitude for Armijo:
    # D = η_alpha * ||grad_alpha||^2 + η_beta * ||grad_beta_free||^2
    D <- eta_alpha * sum(dlda^2) + eta_beta * sum(grad_beta_free^2)
    if (D <= 0) {
      final_iter <- iter - 1L
      break  # gradients essentially zero
    }
    
    # Backtracking line search
    step_scale <- 1.0
    ll_new <- ll_old
    accepted <- FALSE
    
    for (ls_iter in 1:ls_max) {
      
      # Candidate alpha
      alpha_new <- alpha + step_scale * dir_alpha
      
      # Candidate beta free coords
      beta_free_new <- beta[2:Km] + step_scale * dir_beta_free
      
      # Re-enforce beta constraint
      beta_new <- beta
      beta_new[2:Km] <- beta_free_new
      beta_new[1]    <- g_beta(beta_new[2:Km])
      
      ll_candidate <- loglik_fun(alpha_new, beta_new)
      
      # Armijo condition for maximization:
      # ℓ(new) >= ℓ(old) + c * step_scale * D
      if (ll_candidate >= ll_old + ls_c * step_scale * D) {
        ll_new   <- ll_candidate
        alpha    <- alpha_new
        beta     <- beta_new
        accepted <- TRUE
        break
      } else {
        step_scale <- step_scale * ls_rho
      }
    }
    
    # If no step accepted, treat as convergence/stall
    if (!accepted) {
      final_iter <- iter - 1L
      break
    }
    
    # record parameters *after* this accepted update
    param_history[[length(param_history) + 1]] <- rbind(
      data.frame(iter = iter,
                 type = "alpha",
                 k = seq_len(Km),
                 value = alpha),
      data.frame(iter = iter,
                 type = "beta",
                 k = seq_len(Km),
                 value = beta)
    )
    
    # Monotone increase in log-likelihood at accepted steps
    if (abs(ll_new - ll_old) < tol) {
      ll_old    <- ll_new
      final_iter <- iter
      break
    }
    
    ll_old <- ll_new
    final_iter <- iter
  }
  
  # Bind history into a single data.frame
  param_df <- do.call(rbind, param_history)
  
  list(
    alpha   = alpha,
    beta    = beta,
    logLik  = ll_old,
    iter    = final_iter,
    history = param_df
  )
}



out_test <- my_gd_ls(n0=Y[which(X[,2] == 0), ] %>% colSums, 
                     n1 = Y[which(X[,2] == 1), ] %>% colSums, 
                     g_beta=pseudohuber_median_mod,  
                     g_beta_grad= dpseudohuber_median_dx_mod, 
                     eta_alpha = 1e-3,
                     eta_beta  = 1e-3,
                     maxit = 1000,
                     tol = 1e-6)
out_test$alpha
out_test$beta
fn3_1_J$B # j_constr=1, j_ref=J

n0s <- Y[which(X[,2] == 0), ] %>% colSums
n1s <- Y[which(X[,2] == 1), ] %>% colSums
n0s_mod <- n0s; n0s_mod[1] <- n0s[1] + n1s[1]


Y
logit <- function(x) log(x/(1-x))
lg <- logit(n0s/sum(n0s))
lg
lg - lg[J]

logit(n0s_mod/sum(n0s_mod)) - logit(n0s_mod/sum(n0s_mod))[J]


pseudohuber_median_mod(out_test$alpha[-1])
pseudohuber_median_mod(out_test$beta[-1])



out_test <- my_gd_ls(n0=Y[which(X[,2] == 0), ] %>% colSums, 
                     n1 = Y[which(X[,2] == 1), ] %>% colSums, 
                     g_alpha=pseudohuber_median_mod, 
                     g_alpha_grad=dpseudohuber_median_dx_mod,
                     g_beta=pseudohuber_median_mod,  
                     g_beta_grad= dpseudohuber_median_dx_mod, 
                     eta_alpha = 1e-3,
                     eta_beta  = 1e-3,
                     maxit = 1000,
                     tol = 1e-6)
out_test
fn3_1_J$B

pseudohuber_median_mod(out_test$alpha[-1])
pseudohuber_median_mod(out_test$beta[-1])




fit_null(B=fit3$B, Y=Y, X = X, k_constr=2, j_constr=1, j_ref=J,
         constraint_fn=list(pseudohuber_median, pseudohuber_median),
         constraint_grad_fn=list(radEmu::dpseudohuber_median_dx, radEmu::dpseudohuber_median_dx),
         B_tol=1e-8, constraint_tol=1e-8)

fn3_1_J <- fit_null(B=fit3$B, Y=Y, X = X, k_constr=2, j_constr=1, j_ref=J,
                    constraint_fn=list(pseudohuber_median, pseudohuber_median),
                    constraint_grad_fn=list(radEmu::dpseudohuber_median_dx, radEmu::dpseudohuber_median_dx),
                    B_tol=1e-8, constraint_tol=1e-8)
# fn3_1_J$B
# 
# fn3_1_J$B[fn3_1_J$k_constr, setdiff(1:J, fn3_1_J$j_constr)] %>% psuedohuber_median
# fn3_1_J$B[fn3_1_J$k_constr, fn3_1_J$j_constr] 
# 
# fn3_1_J$B[1, ] %>% psuedohuber_median
# fn3_1_J$B[1, ] %>% dput
# fn3_1_J$B[1, setdiff(1:J, c(1, J))] %>% psuedohuber_median
# fn3_1_J$B[1, setdiff(1:J, J)] %>% psuedohuber_median
# fn3_1_J$B[1, setdiff(1:J, 1)] %>% psuedohuber_median
# fn3_1_J$B[1, c(1, J)] %>% psuedohuber_median
# 
# fn3_1_J$B[1, 1:5] %>% psuedohuber_median
# fn3_1_J$B[1, c(2,3,4)] %>% psuedohuber_median
# 
# 
# fn3_1_J$B[1, c(1,5)] %>% psuedohuber_median
# 
# 
# fn3_1_J$B[1, fn3_1_J$j_constr] 
# fn3_1_J$B
# 
# fn3_1_J$B[1, ] %>% psuedohuber_median
