null_repar_one <- function(B,j,z,Y,X,j_constr,k_constr, constraint_fn,constraint_grad_fn,
                           line_search = "armijo",
                           method = "fisher",
                           starting_stepsize = 0.5){
  p <- ncol(X)
  opt_fn <- function(x){ -l_aug(Bj = x[1:p],
                                Bj_constr_no_k_constr = x[(1:(p - 1)) + p],
                                z = z,
                                Y = Y,
                                X = X,
                                Bk_constr_no_j_j_constr = B[k_constr,-c(j,j_constr)],
                                k_constr = k_constr,
                                j = j,
                                j_constr = j_constr,
                                constraint_fn = constraint_fn)}
  
  gr_fn <- function(x){ -d_l_aug_dB(Bj = x[1:p],
                                    Bj_constr_no_k_constr = x[(1:(p - 1)) + p],
                                    z = z,
                                    Y = Y,
                                    X = X,
                                    Bk_constr_no_j_j_constr = B[k_constr,-c(j,j_constr)],
                                    k_constr = k_constr,
                                    j = j,
                                    j_constr = j_constr,
                                    constraint_fn = constraint_fn,
                                    constraint_grad_fn = constraint_grad_fn)}
  
  if(method == "fisher"){
    log_means <- X%*%B[,c(j,j_constr)] + matrix(z,ncol = 1)%*%matrix(1,nrow = 1, ncol = 2)
    
    lil_info <- lapply(1:2,
                       function(j){
                         t(X)%*%diag(exp(log_means[,j]))%*%X
                       })
    lil_info <- Matrix::bdiag(lil_info)
    cg <- constraint_grad_fn(c(B[k_constr,-c(j_constr)],B[k_constr,j]))
    
    pop_index <- p + k_constr
    info_out <- lil_info[pop_index,pop_index]
    lil_info <- lil_info[-pop_index,-pop_index]
    lil_info[k_constr,  k_constr ] <- 
      lil_info[k_constr,  k_constr ] + cg[length(cg)]^2*info_out
    x <- c(B[,j],B[-k_constr,j_constr])
    gr <- gr_fn(x)
    update <- try(-qr.solve(lil_info,gr,tol = 1e-20))
    if(inherits(update,"try-error")){
      update <- -gr
      message("using gradient")
    }
    
    
    
    if(line_search == "exact"){
      message("Using exact line search; j = ",j)
      stepsize <- optim(0,function(eps) opt_fn(x + eps*update),
                        method = "Brent",
                        lower = -starting_stepsize,
                        upper= starting_stepsize)$par
    } else{
      curr_val <- opt_fn(x)
      stepsize <- starting_stepsize
      for_arm <- sum(update*gr)
      
      
      #
      stepsize_ok <- opt_fn(x + stepsize*update) <= curr_val + (1e-4)*stepsize*for_arm
      
      while(!stepsize_ok){
        stepsize <- stepsize*0.25
        stepsize_ok <- opt_fn(x + stepsize*update) <= curr_val + (1e-4)*stepsize*for_arm
        
        if(stepsize < 0.001){
          stepsize_ok <- TRUE
          stepsize <- optim(0,function(eps) opt_fn(x + eps*update),
                            method = "Brent",
                            lower = -starting_stepsize,
                            upper= starting_stepsize)$par
        }
        #
        
      }
    }
  } else{
    
    update <- optim(c(B[,j],B[-k_constr,j_constr]),
                    fn = opt_fn,
                    gr = gr_fn,
                    method = "BFGS",
                    control = list(maxit = 1
                    ))$par - c(B[,j],B[-k_constr,j_constr])
    stepsize <- starting_stepsize
  }

  return(stepsize*update)
}
