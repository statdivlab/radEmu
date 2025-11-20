## this is equivalent to fit_null with k_constr=2, j_constr=1, j_ref=J,
my_fs_stable_two_groups <- function(n0, 
                                    n1,
                                    g_beta,  
                                    g_beta_grad,
                                    maxit = 1000,
                                    tol = 1e-8,
                                    ls_max = 20,       # max step-halving
                                    ls_rho = 0.5) {    # shrink factor for step-halving
  
  K  <- length(n0)
  Km <- K - 1  # m
  
  if (K < 3) {
    stop("beta constraint beta1 = g(beta2,...,beta_{K-1}) requires K >= 3.")
  }
  
  # Parameters: categories 1..K-1 (baseline K is 0)
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
    
    eps <- 1e-12
    p0_full <- pmin(pmax(p0_full, eps), 1 - eps)
    p1_full <- pmin(pmax(p1_full, eps), 1 - eps)
    
    sum(n0 * log(p0_full)) + sum(n1 * log(p1_full))
  }
  
  # ---- history tracking ----
  param_history <- list()
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
  
  ll_old <- loglik_fun(alpha, beta)
  final_iter <- 0L
  
  # tuning constants for stability
  ridge_base     <- 1e-4   # base ridge
  max_step_norm  <- 5      # max Euclidean norm of step in θ
  clip_logit_max <- 15     # clip for alpha/beta components
  
  for (iter in 1:maxit) {
    
    # Probabilities for non-baseline categories 1..K-1
    p0 <- p0fun(alpha)
    p1 <- p1fun(alpha, beta)
    
    # ----- SCORES (gradients) -----
    # Unconstrained score wrt alpha[1:Km] and beta[1:Km]
    s_alpha_full <- (n0[1:Km] + n1[1:Km]) - n0p * p0 - n1p * p1
    s_beta_full  <-  n1[1:Km]            - n1p * p1
    
    # Beta free gradient via chain rule (k = 2..Km)
    grad_g_beta <- g_beta_grad(beta[2:Km])          # length Km-1
    s_beta_free <- s_beta_full[2:Km] + s_beta_full[1] * grad_g_beta
    
    # Score vector for reduced parameters θ = (alpha[1:Km], beta_free[1:(Km-1)])
    score_theta <- c(s_alpha_full, s_beta_free)
    
    # Convergence: score norm
    score_norm <- sqrt(sum(score_theta^2))
    if (score_norm < tol) {
      final_iter <- iter - 1L
      break
    }
    
    # ----- FULL FISHER INFORMATION (α, β_full) -----
    # M0, M1 (m x m)
    M0 <- diag(p0) - tcrossprod(p0)
    M1 <- diag(p1) - tcrossprod(p1)
    
    Iaa <- n0p * M0 + n1p * M1         # alpha-alpha
    Iab <- n1p * M1                    # alpha-beta
    Ibb <- n1p * M1                    # beta-beta
    
    # Full J (2m x 2m)
    J_top    <- cbind(Iaa, Iab)
    J_bottom <- cbind(t(Iab), Ibb)
    J <- rbind(J_top, J_bottom)
    
    # ----- TRANSFORMATION T from θ to φ = (alpha, beta_full) -----
    m  <- Km
    p  <- 2 * m - 1  # dim(θ)
    q  <- 2 * m      # dim(φ)
    
    T <- matrix(0, nrow = q, ncol = p)
    # alpha-part: identity
    T[1:m, 1:m] <- diag(m)
    
    # beta-part: B (m x (m-1)) at rows (m+1)..(2m), cols (m+1)..(2m-1)
    B <- matrix(0, nrow = m, ncol = m - 1)
    B[1, ] <- grad_g_beta                 # row for beta1
    if (m > 1) {
      B[2:m, ] <- diag(m - 1)            # rows beta2..betam
    }
    T[(m + 1):(2 * m), (m + 1):p] <- B
    
    # Fisher info for θ
    I_theta <- crossprod(T, J %*% T)     # T^T J T
    
    # ----- solve for step with adaptive ridge -----
    ridge_multipliers <- c(1, 10, 100, 1000, 1e4)
    step_theta <- rep(NA_real_, length(score_theta))
    solved <- FALSE
    
    for (rm in ridge_multipliers) {
      ridge <- ridge_base * rm
      I_theta_r <- I_theta + ridge * diag(p)
      step_try <- tryCatch(
        solve(I_theta_r, score_theta),
        error = function(e) NA_real_
      )
      if (all(is.finite(step_try))) {
        step_theta <- step_try
        solved <- TRUE
        break
      }
    }
    
    if (!solved) {
      final_iter <- iter - 1L
      break
    }
    
    # Cap step length
    step_norm <- sqrt(sum(step_theta^2))
    if (step_norm > max_step_norm) {
      step_theta <- step_theta * (max_step_norm / step_norm)
    }
    
    # If step is essentially zero, we're done
    if (sqrt(sum(step_theta^2)) < tol) {
      final_iter <- iter - 1L
      break
    }
    
    # ----- Step-halving line search for robustness -----
    step_scale <- 1.0
    ll_new <- ll_old
    accepted <- FALSE
    
    for (ls_iter in 1:ls_max) {
      
      theta_old <- c(alpha, beta[2:Km])
      theta_new <- theta_old + step_scale * step_theta
      
      alpha_new      <- theta_new[1:Km]
      beta_free_new  <- theta_new[(Km + 1):(2 * Km - 1)]
      
      beta_new <- beta
      beta_new[2:Km] <- beta_free_new
      beta_new[1]    <- g_beta(beta_new[2:Km])
      
      # clip logits to avoid numeric blowups
      alpha_new <- pmax(pmin(alpha_new,  clip_logit_max), -clip_logit_max)
      beta_new  <- pmax(pmin(beta_new,   clip_logit_max), -clip_logit_max)
      
      ll_candidate <- loglik_fun(alpha_new, beta_new)
      
      if (is.finite(ll_candidate) && ll_candidate >= ll_old) {
        alpha    <- alpha_new
        beta     <- beta_new
        ll_new   <- ll_candidate
        accepted <- TRUE
        break
      } else {
        step_scale <- step_scale * ls_rho
      }
    }
    
    if (!accepted) {
      final_iter <- iter - 1L
      break
    }
    
    # record parameters after accepted update
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
    
    if (abs(ll_new - ll_old) < tol) {
      ll_old    <- ll_new
      final_iter <- iter
      break
    }
    
    ll_old <- ll_new
    final_iter <- iter
  }
  
  param_df <- do.call(rbind, param_history)
  
  list(
    alpha   = alpha,
    beta    = beta,
    logLik  = ll_old,
    iter    = final_iter,
    history = param_df
  )
}
