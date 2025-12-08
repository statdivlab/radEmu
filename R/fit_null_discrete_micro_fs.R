## this is equivalent to fit_null with j_ref=J
fit_null_discrete_micro_fs <- function(Y, X,
                                       k_constr, j_constr,          # indices: row k* in 1..(p-1), column j* in 1..(J-1)
                                       constraint_fn, constraint_grad_fn,       # g: R^{m-1}->R on the constrained row; grad: R^{m-1}->R^{m-1}
                                       maxit = 1000, 
                                       tol = 1e-8,
                                       ls_max = 20, 
                                       ls_rho = 0.5,
                                       ridge_base = 1e-4,
                                       max_step_norm = 5,
                                       clip_logit_max = 15) {
  Y <- as.matrix(Y); X <- as.matrix(X)
  p <- nrow(Y); J <- ncol(Y); m <- J - 1L
  if (nrow(X) != p || ncol(X) != p - 1L) stop("X must be p x (p-1).")
  if (J < 3) stop("Need J >= 3 (one baseline).")
  if (k_constr < 1L || k_constr > (p - 1L)) stop("k_constr must be in 1..(p-1).")
  if (j_constr < 1L || j_constr > m) stop("j_constr must be in 1..(J-1).")
  
  n <- rowSums(Y)
  
  # Parameters
  alpha <- rep(0, m)
  beta  <- matrix(0, nrow = p - 1L, ncol = m)
  # initialize constrained element coherently (assume g excludes j_constr)
  free_idx <- setdiff(seq_len(m), j_constr)
  beta[k_constr, j_constr] <- constraint_fn(beta[k_constr, free_idx])
  
  softmax_nb <- function(eta) { e <- exp(eta); e / (1 + sum(e)) }
  
  # Log-likelihood
  loglik_fun <- function(alpha, beta) {
    ll <- 0; eps <- 1e-12
    for (i in 1:p) {
      eta <- alpha + drop(X[i, , drop = FALSE] %*% beta) # length m
      pnb <- softmax_nb(eta)
      pfull <- pmin(pmax(c(pnb, 1 - sum(pnb)), eps), 1 - eps)
      ll <- ll + sum(Y[i, ] * log(pfull))
    }
    ll
  }
  
  # History
  param_hist <- list()
  record <- function(iter) {
    rbind(
      data.frame(iter = iter, type = "alpha", row = NA_integer_, col = seq_len(m), value = alpha),
      data.frame(iter = iter, type = "beta",  row = rep(seq_len(p - 1L), each = m),
                 col = rep(seq_len(m), times = p - 1L), value = as.vector(beta))
    )
  }
  param_hist[[1]] <- record(0L)
  
  ll_old <- loglik_fun(alpha, beta); final_iter <- 0L
  
  for (iter in 1:maxit) {
    # probs and M_i
    p_list <- vector("list", p)
    M_list <- vector("list", p)
    for (i in 1:p) {
      eta <- alpha + drop(X[i, , drop = FALSE] %*% beta)
      pnb <- softmax_nb(eta)
      p_list[[i]] <- pnb
      M_list[[i]] <- diag(pnb) - tcrossprod(pnb)
    }
    
    # scores
    s_alpha <- rep(0, m)
    s_beta  <- matrix(0, nrow = p - 1L, ncol = m)
    for (i in 1:p) {
      s_i <- Y[i, 1:m] - n[i] * p_list[[i]] # length m
      s_alpha <- s_alpha + s_i
      # add outer product X[i,]^T * s_i
      s_beta <- s_beta + tcrossprod(X[i, ], s_i) # (p-1) x m
    }
    
    # reduced score: theta = (alpha, rows k != k*, full m; row k* free indices only)
    s_theta <- c(
      s_alpha,
      as.vector(t(s_beta[setdiff(seq_len(p - 1L), k_constr), , drop = FALSE])),
      # constrained row (chain rule)
      {
        grad_g <- constraint_grad_fn(beta[k_constr, free_idx])
        s_beta[k_constr, free_idx] + s_beta[k_constr, j_constr] * grad_g
      }
    )
    
    if (sqrt(sum(s_theta^2)) < tol) { final_iter <- iter - 1L; break }
    
    # Fisher J in full phi = (alpha, vec(beta rows 1..p-1 each length m))
    q_dim <- m + (p - 1L) * m
    Jmat <- matrix(0, q_dim, q_dim)
    
    # indices
    alpha_idx <- 1:m
    row_block_idx <- function(k) { # k in 1..(p-1), returns indices for row k in vec(beta)
      start <- m + (k - 1L) * m + 1L
      start:(start + m - 1L)
    }
    
    for (i in 1:p) {
      Mi <- M_list[[i]]; ni <- n[i]; Xi <- X[i, , drop = FALSE]          # 1 x (p-1)
      # alpha-alpha
      Jmat[alpha_idx, alpha_idx] <- Jmat[alpha_idx, alpha_idx] + ni * Mi
      # alpha-beta rows
      for (k in 1:(p - 1L)) {
        rb <- row_block_idx(k)
        Jmat[alpha_idx, rb] <- Jmat[alpha_idx, rb] + ni * (Mi * Xi[1, k])
        Jmat[rb, alpha_idx] <- t(Jmat[alpha_idx, rb])
      }
      # beta-beta rows
      for (k in 1:(p - 1L)) {
        for (l in 1:(p - 1L)) {
          rbk <- row_block_idx(k); rbl <- row_block_idx(l)
          Jmat[rbk, rbl] <- Jmat[rbk, rbl] + ni * (Mi * (Xi[1, k] * Xi[1, l]))
        }
      }
    }
    
    # Transformation T: reduced theta -> full phi
    # theta packs: alpha, then rows k != k* (full m), then row k* (m-1 free entries)
    p_dim <- m + (p - 2L) * m + (m - 1L)
    T <- matrix(0, nrow = q_dim, ncol = p_dim)
    # alpha identity
    T[alpha_idx, 1:m] <- diag(m)
    col_ptr <- m + 1L
    # rows k != k*
    for (k in setdiff(seq_len(p - 1L), k_constr)) {
      rb <- row_block_idx(k)
      T[rb, col_ptr:(col_ptr + m - 1L)] <- diag(m)
      col_ptr <- col_ptr + m
    }
    # row k* mapping B (m x (m-1))
    grad_g <- constraint_grad_fn(beta[k_constr, free_idx])
    Bmap <- matrix(0, nrow = m, ncol = m - 1L)
    # identity on free indices
    for (pos in seq_along(free_idx)) {
      Bmap[free_idx[pos], pos] <- 1
    }
    # constrained element row = grad_g^T
    Bmap[j_constr, ] <- grad_g
    T[row_block_idx(k_constr), col_ptr:(col_ptr + (m - 2L))] <- Bmap
    
    I_theta <- crossprod(T, Jmat %*% T)
    
    # Adaptive ridge solve
    ridge_mult <- c(1, 10, 100, 1000, 1e4)
    step_theta <- rep(NA_real_, p_dim); solved <- FALSE
    for (rm in ridge_mult) {
      Ir <- I_theta + (ridge_base * rm) * diag(p_dim)
      st <- tryCatch(solve(Ir, s_theta), error = function(e) NA_real_)
      if (all(is.finite(st))) { step_theta <- st; solved <- TRUE; break }
    }
    if (!solved) { final_iter <- iter - 1L; break }
    
    # Step capping
    sn <- sqrt(sum(step_theta^2))
    if (sn > max_step_norm) step_theta <- step_theta * (max_step_norm / sn)
    if (sqrt(sum(step_theta^2)) < tol) { final_iter <- iter - 1L; break }
    
    # Pack current reduced theta
    pack_theta <- function(alpha, beta) {
      c(
        alpha,
        as.vector(t(beta[setdiff(seq_len(p - 1L), k_constr), , drop = FALSE])),
        beta[k_constr, free_idx]
      )
    }
    unpack_theta <- function(vec) {
      wp <- 1L
      a_new <- vec[wp:(wp + m - 1L)]; wp <- wp + m
      beta_new <- beta
      # rows k != k*
      for (k in setdiff(seq_len(p - 1L), k_constr)) {
        beta_new[k, ] <- vec[wp:(wp + m - 1L)]; wp <- wp + m
      }
      # row k* free
      bfree <- vec[wp:(wp + (m - 2L))]
      beta_new[k_constr, free_idx] <- bfree
      beta_new[k_constr, j_constr]   <- constraint_fn(bfree)
      list(alpha = a_new, beta = beta_new)
    }
    
    # Line search (monotone ℓ)
    theta_old <- pack_theta(alpha, beta)
    step_scale <- 1.0; accepted <- FALSE; ll_new <- ll_old
    for (ls in 1:ls_max) {
      prop <- unpack_theta(theta_old + step_scale * step_theta)
      alpha_new <- prop$alpha
      beta_new  <- prop$beta
      # clip logits (keep numerics sane)
      alpha_new <- pmax(pmin(alpha_new,  clip_logit_max), -clip_logit_max)
      beta_new  <- pmax(pmin(beta_new,   clip_logit_max), -clip_logit_max)
      ll_cand <- loglik_fun(alpha_new, beta_new)
      if (is.finite(ll_cand) && ll_cand >= ll_old) {
        alpha <- alpha_new; beta <- beta_new; ll_new <- ll_cand; accepted <- TRUE; break
      }
      step_scale <- step_scale * ls_rho
    }
    if (!accepted) { final_iter <- iter - 1L; break }
    
    # record and check convergence
    param_hist[[length(param_hist) + 1L]] <- record(iter)
    if (abs(ll_new - ll_old) < tol) { ll_old <- ll_new; final_iter <- iter; break }
    ll_old <- ll_new; final_iter <- iter
  }
  
  B <- cbind(rbind(alpha, beta), 0)
  
  history <- do.call(rbind, param_hist)
  list(B = B,                 
       logLik = ll_old,
       iter   = final_iter,
       history = history)
}
