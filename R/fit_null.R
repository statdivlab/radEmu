fit_null <- function(B,
                     Y, #Y (with augmentations)
                     X, #design matrix
                     X_cup = NULL,
                     k_constr, #row index of B to constrain
                     j_constr, #col index of B to constrain
                     j_ref,
                     constraint_fn, #constraint function
                     constraint_grad_fn, #gradient of constraint fn\
                     rho_init = 1,
                     tau = 1.2,
                     kappa = 0.8,
                     B_tol = 1e-2,
                     inner_tol = 0.01,
                     constraint_tol = 1e-4,
                     max_step = 5,
                     c1 = 1e-4,
                     maxit = 1000,
                     inner_maxit = 25,
                     n_final_inner_it = 25,
                     verbose = FALSE,
                     trackB = FALSE
) {
  
  if (trackB) warning("trackB needs fixing; pieces are missing")
  
  J <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  
  #store cols of B to update in a vector
  j_to_update <- (1:J)[(1:J) != j_ref]
  #in vector format, which indexes to remove (bc B^j_constr = 0 id. constr.)
  indexes_to_remove <- (j_ref - 1)*p + 1:p
  
  
  #change id. constr. by subtracting approp. col of B from others
  for(k in 1:p) {
    B[k,] <- B[k,] - B[k,j_ref]
  }
  
  #how large is norm of score at initiation?
  z <- update_z(X = X, Y = Y, B = B)
  
  log_mean <- X%*%B +
    matrix(z,ncol = 1)%*%matrix(1,ncol = J, nrow = 1)
  
  # if (trackB) {
  #   Bs <- data.frame(k = rep(1:p,each = J),
  #                    j = rep(1:J,p),
  #                    B = NA)
  # } else{
  Bs <- NULL
  # }
  
  #get X_cup for later use
  if (is.null(X_cup)) {
    X_cup = X_cup_from_X(X,J)
  }
  
  #set  <- ation to zero
  iter <- 0
  
  #compute gap
  gap <- constraint_fn(B[k_constr,]) - B[k_constr,j_constr]
  init_gap <- abs(gap)
  
  #set rho equal to initial value
  rho <- rho_init
  
  #initiate u
  u <- rho*gap
  
  # if (trackB) {
  #   Bs$iter <- 0
  #   Bs$inner_iter <- 0
  #   Bs$rho <- rho
  #   Bs$u <- u
  #   Bs$Q <- NA
  #   Bs$gap <- gap
  #   Bs$score_stat <- 0
  #   Bs$score_err <- 0
  #   Bs$dlag_dB <- NA
  #   lilBs <- Bs
  #   
  #   for(k in 1:p) {
  #     Bs$B[1:J + J*(k - 1)] <- B[k,]
  #   }
  # }
  
  loop_j_to_update <- j_to_update
  
  use_max_inner_it <- FALSE
  B_diff <- Inf
  
  while((abs(gap) > constraint_tol | B_diff> B_tol #proj_score > B_tol #| Q > inner_tol
  ) & iter <= maxit) {
    
    iter <- iter + 1    #increment iteration
    inner_iter <- 0     #initiate internal iteration
    
    #evaluate augmented Lagrangian
    log_means <- do.call(cbind,lapply(1:J,function(j) X%*%B[,j] + z))
    
    #     message("Have to update how lagrangian and its derivs are computed to reflect
    # better ref taxon")
    curr_lag_val <- sum(-log_means*Y + exp(log_means)) + u*gap + (rho/2)*(gap^2)
    
    old_gap <- gap
    old_B <- B
    inner_diff <- Inf
    while(((inner_diff > inner_tol | use_max_inner_it) & inner_iter <= inner_maxit)| inner_iter == 0) {
      
      inner_old_B <- B
      
      update <- macro_fisher_null(X = X,
                                  Y = Y,
                                  B = B,
                                  z = z,
                                  J = J,
                                  p = p,
                                  k_constr = k_constr,
                                  j_constr = j_constr,
                                  j_ref = j_ref,
                                  rho = rho,
                                  u = u,
                                  max_step = max_step,
                                  constraint_fn = constraint_fn,
                                  constraint_grad_fn = constraint_grad_fn,
                                  verbose = verbose, 
                                  c1 = 1e-4)
      B <- B + update$update
      gap <- update$gap
      z <- update_z(Y,X,B)
      
      # if (trackB) {
      #   
      #   ### AW notes 1/15/24 that dlag_lamhat_dB is not computed in this function. 
      #   
      #   temp_Bs <- lilBs
      #   for(k in 1:p) {
      #     temp_Bs$B[1:J + J*(k - 1)] <- B[k,]
      #   }
      #   temp_Bs$iter <- iter
      #   temp_Bs$inner_iter <- inner_iter
      #   temp_Bs$rho <- rho
      #   temp_Bs$u <- u
      #   temp_Bs$Q <- NA
      #   temp_Bs$gap <- gap
      #   if (iter > 1 | inner_iter > 1) {
      #     jcounter <- 1
      #     for(j in 1:J) {
      #       if (j == j_constr) {
      #         temp_Bs$dlag_dB[temp_Bs$j ==j] <- 0
      #       } else{
      #         
      #         ### maybe in 
      #         
      #         # dlag_lamhat_dB <- dldB + lambda_hat*dgdB
      #         # temp_Bs$dlag_dB[temp_Bs$j ==j] <- dlag_lamhat_dB[(jcounter - 1)*p + 1:p,]
      #         # jcounter <- jcounter + 1
      #       }
      #     }
      #   }
      #   
      #   Bs <- rbind(Bs,temp_Bs)
      #   
      # }
      
      inner_diff <- max(abs(B - inner_old_B))
      
      inner_iter <- inner_iter + 1
      
    }
    
    B_diff <- max(abs(B - old_B))
    
    if (verbose) {
      message("Max absolute difference in B since last augmented Lagrangian outer step: ", signif (B_diff,3))
    }
    
    
    if ( abs(gap) > constraint_tol) {
      u <- u + rho*gap
      if (abs(gap/old_gap) > kappa) {
        rho <- rho*tau
      } 
    }
    
    if (verbose) {
      message("Estimate deviates from feasibility by ",signif (gap,3))
      message("Parameter u set to ", signif (u,3))
      message("Parameter rho set to ", signif (rho,3))
    }
    
    if (abs(gap) < constraint_tol) {
      use_max_inner_it <- TRUE
      # message("Score stat: ",signif (check$score_stat*(n/(n - 1)),3))
      # message("Influence stat: ",signif (check$influence_stat,3))
      # proj_score <- as.numeric(as.matrix(check$influence_stat))
    } else{
      use_max_inner_it <- FALSE
    }
    
  }
  
  return(list("B" = B,
              "k_constr" = k_constr,
              "j_constr" = j_constr,
              # "score_stats" = list(score_stats),
              # "proj_score" = proj_score,
              "niter" = iter,
              "gap" = gap,
              "u" = u,
              "rho" = rho,
              "Bs" = Bs))
}

