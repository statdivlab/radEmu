my_gd_fs <- function(n0, n1,
                     g_beta,  g_beta_grad,
                     eta_alpha = 1e-3,  # kept for signature; not used
                     eta_beta  = 1e-3,  # kept for signature; not used
                     maxit = 1000,
                     B_tol = 1e-8,
                     ls_max = 20,       # step-halving max
                     ls_rho = 0.5) {    # step shrink factor
  
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
  
  for (iter in 1:maxit) {
    
    # Probabilities for non-baseline categories 1..K-1
    p0 <- p0fun(alpha)
    p1 <- p1fun(alpha, beta)
    
    # ----- SCORES (gradients) -----
    # Unconstrained score wrt alpha[1:Km] and beta[1:Km]
    s_alpha_full <- (n0[1:Km] + n1[1:Km]) - n0p * p0 - n1p * p1
    s_beta_full  <-  n1[1:Km]            - n1p * p1
    
    # Beta free gradient via chain rule (k = 2..Km)
    grad_g_beta     <- g_beta_grad(beta[2:Km])          # length Km-1
    s_beta_free     <- s_beta_full[2:Km] + s_beta_full[1] * grad_g_beta
    
    # Score vector for reduced parameters θ = (alpha[1:Km], beta_free[1:(Km-1)])
    score_theta <- c(s_alpha_full, s_beta_free)
    
    # Convergence: score norm
    score_norm <- sqrt(sum(score_theta^2))
    # if (score_norm < B_tol) {
    #   final_iter <- iter - 1L
    #   break
    # }
    
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
    
    # Add tiny ridge for numerical stability if needed
    ridge <- 1e-8
    I_theta_r <- I_theta + ridge * diag(p)
    
    # Fisher scoring step: θ_new = θ_old + I^{-1} S
    step_theta <- tryCatch(
      solve(I_theta_r, score_theta),
      error = function(e) {
        # In case of numerical issues, break
        rep(0, length(score_theta))
      }
    )
    
    # If step is essentially zero, we're done
    # if (sqrt(sum(step_theta^2)) < B_tol) {
    if (max(abs(step_theta)) < B_tol) {
      final_iter <- iter - 1L
      break
    }
    
    # ----- Step-halving line search for robustness -----
    step_scale <- 1.0
    ll_new <- ll_old
    accepted <- FALSE
    
    for (ls_iter in 1:ls_max) {
      
      theta_new <- c(alpha, beta[2:Km]) + step_scale * step_theta
      alpha_new <- theta_new[1:Km]
      beta_free_new <- theta_new[(Km + 1):(2 * Km - 1)]
      
      beta_new <- beta
      beta_new[2:Km] <- beta_free_new
      beta_new[1]    <- g_beta(beta_new[2:Km])
      
      print(cbind(alpha_new, beta_new))
      ll_candidate <- loglik_fun(alpha_new, beta_new)
      print(ll_candidate)
      
      if (ll_candidate >= ll_old) {
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
    
    if (abs(ll_new - ll_old) < B_tol) {
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
