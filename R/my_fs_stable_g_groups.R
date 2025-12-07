my_fs_stable_g_groups <- function(n_list,
                                  g_beta,  
                                  g_beta_grad,
                                  maxit = 1000,
                                  tol = 1e-8,
                                  ls_max = 20,       # max step-halving
                                  ls_rho = 0.5) {    # shrink factor for step-halving
  
  G <- length(n_list)        # number of groups
  if (G < 2) stop("Need at least 2 groups (baseline + at least 1 non-baseline).")
  
  K <- length(n_list[[1]])
  if (K < 3) stop("beta1 = g(beta2,...,beta_{K-1}) requires K >= 3.")
  for (g in seq_len(G)) {
    if (length(n_list[[g]]) != K) stop("All n_g must have same length K.")
  }
  
  m <- K - 1  # non-baseline categories
  
  # Group totals
  n_tot <- vapply(n_list, sum, numeric(1))
  
  # Parameters:
  # alpha: length m
  # beta_mat: (G-1) x m, each row = beta^(g) for group g=2..G
  alpha   <- rep(0, m)
  beta_mat <- matrix(0, nrow = G - 1, ncol = m)
  
  # Enforce beta constraint initially: for each non-baseline group
  if (G > 1) {
    for (r in 1:(G - 1)) {
      beta_mat[r, 1] <- g_beta(beta_mat[r, 2:m])
    }
  }
  
  # Softmax for non-baseline categories
  softmax_nb <- function(eta) {
    e <- exp(eta)
    e / (1 + sum(e))
  }
  
  # Log-likelihood given alpha, beta_mat
  loglik_fun <- function(alpha, beta_mat) {
    ll <- 0
    eps <- 1e-12
    
    for (g in seq_len(G)) {
      n_g  <- n_list[[g]]
      if (g == 1) {
        eta <- alpha
      } else {
        r <- g - 1
        eta <- alpha + beta_mat[r, ]
      }
      p_nb <- softmax_nb(eta)
      p_full <- c(p_nb, 1 - sum(p_nb))
      p_full <- pmin(pmax(p_full, eps), 1 - eps)
      ll <- ll + sum(n_g * log(p_full))
    }
    ll
  }
  
  # ---- history tracking ----
  param_history <- list()
  # iteration 0 (initial)
  {
    df_alpha <- data.frame(iter = 0L,
                           type  = "alpha",
                           group = NA_integer_,
                           k     = seq_len(m),
                           value = alpha)
    df_beta_list <- vector("list", G - 1)
    for (r in 1:(G - 1)) {
      df_beta_list[[r]] <- data.frame(iter = 0L,
                                      type  = "beta",
                                      group = r + 1L,  # group index in n_list (2..G)
                                      k     = seq_len(m),
                                      value = beta_mat[r, ])
    }
    param_history[[1]] <- rbind(df_alpha, do.call(rbind, df_beta_list))
  }
  # ---------------------------
  
  ll_old <- loglik_fun(alpha, beta_mat)
  final_iter <- 0L
  
  # stability parameters
  ridge_base     <- 1e-4
  max_step_norm  <- 5
  clip_logit_max <- 15
  
  for (iter in 1:maxit) {
    
    # --- compute p_g and M_g for each group ---
    p_list <- vector("list", G)
    M_list <- vector("list", G)
    for (g in seq_len(G)) {
      if (g == 1) {
        eta <- alpha
      } else {
        r <- g - 1
        eta <- alpha + beta_mat[r, ]
      }
      p_nb <- softmax_nb(eta)
      p_list[[g]] <- p_nb
      M_list[[g]] <- diag(p_nb) - tcrossprod(p_nb)
    }
    
    # --- SCORES ---
    # alpha score
    s_alpha <- numeric(m)
    for (g in seq_len(G)) {
      n_g  <- n_list[[g]]
      n_gp <- n_tot[g]
      p_nb <- p_list[[g]]
      s_alpha <- s_alpha + (n_g[1:m] - n_gp * p_nb)
    }
    
    # beta scores: full and free for each non-baseline group
    s_beta_full_mat <- matrix(0, nrow = G - 1, ncol = m)
    s_beta_free_mat <- matrix(0, nrow = G - 1, ncol = m - 1)
    grad_g_beta_mat <- matrix(0, nrow = G - 1, ncol = m - 1)
    
    if (G > 1) {
      for (r in 1:(G - 1)) {
        g <- r + 1
        n_g  <- n_list[[g]]
        n_gp <- n_tot[g]
        p_nb <- p_list[[g]]
        
        s_beta_full <- n_g[1:m] - n_gp * p_nb
        s_beta_full_mat[r, ] <- s_beta_full
        
        grad_g <- g_beta_grad(beta_mat[r, 2:m])  # length m-1
        grad_g_beta_mat[r, ] <- grad_g
        
        s_beta_free_mat[r, ] <- s_beta_full[2:m] + s_beta_full[1] * grad_g
      }
    }
    
    # Reduced score vector θ = (alpha, all beta_free)
    score_theta <- c(s_alpha, as.vector(t(s_beta_free_mat)))
    
    score_norm <- sqrt(sum(score_theta^2))
    if (score_norm < tol) {
      final_iter <- iter - 1L
      break
    }
    
    # --- FULL FISHER J for (alpha, all beta_full) ---
    q_dim <- m + (G - 1) * m
    J <- matrix(0, nrow = q_dim, ncol = q_dim)
    
    # indices for alpha in full φ
    alpha_idx <- 1:m
    
    # helper to get beta block indices for group r (r=1..G-1)
    beta_block_idx <- function(r) {
      start <- m + (r - 1) * m + 1
      end   <- m + r * m
      start:end
    }
    
    for (g in seq_len(G)) {
      M_g  <- M_list[[g]]
      n_gp <- n_tot[g]
      
      # alpha-alpha block
      J[alpha_idx, alpha_idx] <- J[alpha_idx, alpha_idx] + n_gp * M_g
      
      if (g > 1) {
        r <- g - 1
        B_idx <- beta_block_idx(r)
        
        # alpha-beta and beta-beta blocks
        J[alpha_idx, B_idx] <- J[alpha_idx, B_idx] + n_gp * M_g
        J[B_idx, alpha_idx] <- t(J[alpha_idx, B_idx])
        J[B_idx, B_idx]     <- J[B_idx, B_idx]     + n_gp * M_g
      }
    }
    
    # --- TRANSFORMATION T from θ to φ ---
    p_dim <- m + (G - 1) * (m - 1)
    T <- matrix(0, nrow = q_dim, ncol = p_dim)
    
    # alpha part: identity
    T[alpha_idx, 1:m] <- diag(m)
    
    # beta parts: each group r has a B^(r)
    for (r in 1:(G - 1)) {
      B_idx <- beta_block_idx(r)
      C_start <- m + (r - 1) * (m - 1) + 1
      C_end   <- m + r * (m - 1)
      C_idx   <- C_start:C_end
      
      B <- matrix(0, nrow = m, ncol = m - 1)
      B[1, ] <- grad_g_beta_mat[r, ]
      if (m > 1) {
        B[2:m, ] <- diag(m - 1)
      }
      T[B_idx, C_idx] <- B
    }
    
    I_theta <- crossprod(T, J %*% T)
    
    # --- solve for Fisher scoring step with adaptive ridge ---
    ridge_base <- 1e-4
    ridge_multipliers <- c(1, 10, 100, 1000, 1e4)
    step_theta <- rep(NA_real_, length(score_theta))
    solved <- FALSE
    
    for (rm in ridge_multipliers) {
      ridge <- ridge_base * rm
      I_theta_r <- I_theta + ridge * diag(p_dim)
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
    
    if (sqrt(sum(step_theta^2)) < tol) {
      final_iter <- iter - 1L
      break
    }
    
    # --- Step-halving line search ---
    step_scale <- 1.0
    ll_new <- ll_old
    accepted <- FALSE
    
    for (ls_iter in 1:ls_max) {
      
      theta_old <- c(alpha, as.vector(t(beta_mat[, 2:m])))
      theta_new <- theta_old + step_scale * step_theta
      
      # unpack theta_new
      alpha_new <- theta_new[1:m]
      if (G > 1) {
        bfree_vec <- theta_new[(m + 1):p_dim]
        beta_free_mat_new <- matrix(bfree_vec, nrow = G - 1, ncol = m - 1, byrow = TRUE)
      } else {
        beta_free_mat_new <- matrix(0, nrow = 0, ncol = m - 1)
      }
      
      beta_mat_new <- beta_mat
      if (G > 1) {
        beta_mat_new[, 2:m] <- beta_free_mat_new
        # re-enforce constraint for each non-baseline group
        for (r in 1:(G - 1)) {
          beta_mat_new[r, 1] <- g_beta(beta_mat_new[r, 2:m])
        }
      }
      
      # clip logits to avoid blowups
      alpha_new  <- pmax(pmin(alpha_new,  clip_logit_max), -clip_logit_max)
      beta_mat_new <- pmax(pmin(beta_mat_new, clip_logit_max), -clip_logit_max)
      
      ll_candidate <- loglik_fun(alpha_new, beta_mat_new)
      
      if (is.finite(ll_candidate) && ll_candidate >= ll_old) {
        alpha    <- alpha_new
        beta_mat <- beta_mat_new
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
    df_alpha <- data.frame(iter = iter,
                           type  = "alpha",
                           group = NA_integer_,
                           k     = seq_len(m),
                           value = alpha)
    df_beta_list <- vector("list", G - 1)
    for (r in 1:(G - 1)) {
      df_beta_list[[r]] <- data.frame(iter = iter,
                                      type  = "beta",
                                      group = r + 1L,
                                      k     = seq_len(m),
                                      value = beta_mat[r, ])
    }
    param_history[[length(param_history) + 1]] <- rbind(df_alpha, do.call(rbind, df_beta_list))
    
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
    beta    = beta_mat,
    logLik  = ll_old,
    iter    = final_iter,
    history = param_df
  )
}
