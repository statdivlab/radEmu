
#block approximate Newton update to B used in 
#optimization under null

macro_fisher_null <- function(X,
                              Y,
                              B,
                              z,
                              J,
                              p,
                              k_constr,
                              j_constr,
                              j_ref,
                              rho,
                              u,
                              constraint_fn,
                              constraint_grad_fn,
                              stepsize = 0.5,
                              c1 = 1e-4,
                              regularization = 0.1,
                              max_step = 1,
                              verbose = FALSE,
                              debug = FALSE #return more information useful for
                                            #debugging
                              ){
  #list of log means... by taxon
  log_means_by_j <-
    lapply(1:J,function(j)
      X%*%B[,j,drop = FALSE] + z)
  
  #compute information matrix of *Poisson* log likelihood
  #for fixed z
  info_bits <- lapply(log_means_by_j,
                      function(log_mu){
                        crossprod(X,diag(as.numeric(exp(log_mu))))%*%X
                      })
  info_diag <- do.call(c,lapply(info_bits,diag))
  mean_diag <- mean(info_diag)
  
  info_bits <- lapply(info_bits,
                      function(x){
                        diag(x) <- diag(x) + regularization*mean_diag;
                        return(x)
                      })
  #invert info by inverting each block diagonal element
  inv_info_by_j <- lapply(info_bits,
                          function(x)
                            qr.solve(x,tol = 1e-20))
  #remove diagonal block corresponding to convenience constraint
  inv_info_by_j <- inv_info_by_j[-j_ref]
  
  #stick the block diagonal elements together into one matrix to get 
  #inverse as an actual matrix
  info_inverse <- Matrix::bdiag(inv_info_by_j)
  
  #compute gradient of constraint function
  cg <- constraint_grad_fn(B[k_constr,])
  cg[j_constr] <- cg[j_constr] - 1
  cg <- cg[-j_ref]
  
  #vector of p elements whose k_constr-th element is 1, all else 0
  e_k_constr <- Matrix(0,nrow = p,ncol = 1)
  e_k_constr[k_constr,] <- 1
  #take kronecker producr with gradient of constraint to 
  #get a gradient in long format B (i.e., B_cup)
  cg_expanded <- Matrix::kronecker(Matrix(cg,ncol = 1),e_k_constr)
  v <- sqrt(rho)*cg_expanded
  #do sherman-morrison rank-1 update in a couple steps
  sm_denom <- as.numeric(as.matrix(1 + Matrix::crossprod(v,info_inverse)%*%v))
  sm_half_num <- info_inverse%*%v
  
  #strictly speaking, info_inverse isn't the inverse info -- its the inverse of 
  # (an approximation to) the hessian of the augmented lagrangian
  
  #compute derivative of augmented lagrangian
  lag_deriv <- lapply(1:J,
                      function(j){
                        
                        -t(X)%*%(Y[,j,drop = FALSE] - exp(log_means_by_j[[j]]))
                      })
  lag_deriv <- do.call(rbind,lag_deriv[-j_ref])
  gap <- constraint_fn(B[k_constr,]) - B[k_constr,j_constr]
  lag_deriv <- lag_deriv + (u + rho*gap)*cg_expanded
  
  
  #direction for (approximate) Newton step
  update_dir <- -info_inverse %*%lag_deriv + 
    sm_half_num%*% Matrix::crossprod(sm_half_num,lag_deriv)/sm_denom
  #shorten step if any of its elements are larger than allowed by max_step
  if(max(abs(update_dir))>max_step){
    update_dir <- update_dir/max(abs(update_dir))
  }
  armijo_term <- c1*sum(update_dir*lag_deriv)

  #update direction in long format
  update_dir <- B_from_B_cup(update_dir,J - 1, p)
  #don't update betas for j_ref
  update_dir <- do.call(cbind,list(update_dir[,(1:(J - 1))<j_ref,drop = FALSE],
                                   0,
                                   update_dir[,(1:(J - 1))>=j_ref,drop = FALSE]))
  
  
  #in debug mode, compute full inverse info for backward compatibility
  if(debug){
    info_inverse <- info_inverse - Matrix::tcrossprod(sm_half_num)/sm_denom
  }
  
  
  
  
  
  #find stepsize that satisfies armijo rule
  accept <- FALSE
  
  curr_lag_val <- do.call(sum,
                          lapply(1:J, function(j)
                            -Y[,j,drop = FALSE]*log_means_by_j[[j]] + exp(log_means_by_j[[j]]))) +
    u*gap + (rho/2)*gap^2
  
  while(!accept){
    prop_B <- B + stepsize*update_dir
    prop_gap <- constraint_fn(prop_B[k_constr,]) - prop_B[k_constr,j_constr]
    
    prop_log_means_by_j <-
      lapply(1:J,function(j)
        X%*%prop_B[,j,drop = FALSE] + z)
    
    prop_lag_val <- do.call(sum,
                            lapply(1:J, function(j)
                              -Y[,j,drop = FALSE]*prop_log_means_by_j[[j]] + exp(prop_log_means_by_j[[j]]))) +
      u*prop_gap + (rho/2)*prop_gap^2
    
    accept <- prop_lag_val <= curr_lag_val + armijo_term*stepsize
    stepsize <- stepsize/2
    # stepsize <- stepsize*0.9
    
    # steps <- c(steps,stepsize)
    # props <- c(props,prop_lag_val)
    
  }
  stepsize <- stepsize*2
  if (stepsize <  1e-2 & verbose) {
    #ideally, stepsize doesn't get super small and we never see this message
    #but it can be helpful if something is going wrong
    message("Small stepsize: ",signif(stepsize,2))
  }
  if (debug) {
    return(list(update =  update_dir*stepsize,gap = prop_gap,
                lag_deriv = lag_deriv,
                info_inverse = info_inverse,
                stepsize = stepsize))
  } else{
    return(list(update = update_dir*stepsize,gap = prop_gap))
    
  }
  
  
}
