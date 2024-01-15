

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
                              debug = FALSE){
  log_means_by_j <-
    lapply(1:J,function(j)
      X%*%B[,j,drop = FALSE] + z)
  
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
  
  inv_info_by_j <- lapply(info_bits,
                          function(x)
                            qr.solve(x,tol = 1e-20))
  inv_info_by_j <- inv_info_by_j[-j_ref]
  
  info_inverse <- Matrix::bdiag(inv_info_by_j)
  
  cg <- constraint_grad_fn(B[k_constr,])
  cg[j_constr] <- cg[j_constr] - 1
  cg <- cg[-j_ref]
  
  e_k_constr <- Matrix(0,nrow = p,ncol = 1)
  e_k_constr[k_constr,] <- 1
  cg_expanded <- Matrix::kronecker(Matrix(cg,ncol = 1),e_k_constr)
  v <- sqrt(rho)*cg_expanded
  sm_denom <- as.numeric(as.matrix(1 + Matrix::crossprod(v,info_inverse)%*%v))
  sm_half_num <- info_inverse%*%v
  # sm_half_num <- Matrix(round(sm_half_num,-floor(log(max(abs(sm_half_num)))/log(10) -1)))
  # sm_half_num <- methods::as(sm_half_num,"sparseMatrix")
  info_inverse <- info_inverse - Matrix::tcrossprod(sm_half_num)/sm_denom
  
  # which_addl <- which(cg/sum(cg)>= 1/J)
  # n_addl <- length(which_addl)
  # C <- Matrix(0,ncol = n_addl,nrow = n_addl)
  # for(i in 1:n_addl){
  #   for(j in i:n_addl){
  #     C[i,j] <-
  #       C[j,i] <-
  #       constraint_hess_fn(B[k_constr,],which_addl[i],which_addl[j])
  #   }
  # }
  # U <- Matrix::kronecker(Matrix(as.numeric((1:(J - 1)) %in% which_addl),ncol = 1),e_k_constr)
  # V <- Matrix::t(U)*(u + rho*gap)
  # woodbury_center <- Matrix::solve(Matrix::solve(C) + V%*%info_inverse%*%U)
  #
  # info_inverse <- info_inverse - info_inverse%*%U%*%woodbury_center%*%V%*%info_inverse
  
  ######
  #
  # lag_fn <- function(B_vec_no_ref){
  #   B_no_ref <- B_from_B_cup(B_vec_no_ref,J = J -1, p = p)
  #   B_with_ref <- cbind(B_no_ref[,1:(j_ref - 1)],0,
  #                       B_no_ref[,(j_ref):(J - 1)])
  #   log_means <- lapply(1:J,
  #                       function(x) X%*%B_with_ref[,x,drop = FALSE] + z)
  #   log_means <- do.call(cbind,log_means)
  #
  #   neg_ll <- -sum(Y*log_means - exp(log_means))
  #
  #   fn_gap <- constraint_fn(B_with_ref[k_constr,]) - B_with_ref[k_constr,j_constr]
  #
  #   return(neg_ll + u*fn_gap + (rho/2)*(fn_gap^2))
  # }
  #
  # nderiv <- numDeriv::grad(lag_fn,
  #                          x = as.numeric(as.matrix(B_cup_from_B(B[,-j_ref]))))
  #
  # ninfo <- numDeriv::hessian(lag_fn,
  #                            x = as.numeric(as.matrix(B_cup_from_B(B[,-j_ref]))))
  
  ######
  lag_deriv <- lapply(1:J,
                      function(j){
                        
                        -t(X)%*%(Y[,j,drop = FALSE] - exp(log_means_by_j[[j]]))
                      })
  lag_deriv <- do.call(rbind,lag_deriv[-j_ref])
  gap <- constraint_fn(B[k_constr,]) - B[k_constr,j_constr]
  lag_deriv <- lag_deriv + (u + rho*gap)*cg_expanded
  
  ####
  # plot(asinh(nderiv),asinh(lag_deriv))
  # abline(a = 0, b= 1,lty = 2, col = "red")
  #
  # plot(asinh(diag(ninfo)),
  #      asinh(diag(solve(info_inverse))))
  #
  #
  # abline(a = 0, b= 1,lty = 2, col = "red")
  
  ####
  
  update_dir <- -info_inverse%*%lag_deriv
  # update_dir <- lag_deriv
  if(max(abs(update_dir))>max_step){
    update_dir <- update_dir/max(abs(update_dir))
  }
  armijo_term <- c1*sum(update_dir*lag_deriv)
  # print(armijo_term)
  update_dir <- B_from_B_cup(update_dir,J - 1, p)
  
  update_dir <- do.call(cbind,list(update_dir[,(1:(J - 1))<j_ref,drop = FALSE],
                                   0,
                                   update_dir[,(1:(J - 1))>=j_ref,drop = FALSE]))
  
  
  
  
  
  
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
