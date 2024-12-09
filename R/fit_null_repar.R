#' fits model with B_kj constrained to equal g(B_k) for constraint fn g
#' 
#' @param B description
#' @param Y Y (with augmentations)
#' @param X design matrix
#' @param k_constr row index of B to constrain
#' @param j_constr col index of B to constrain
#' @param j_ref column index of convenience constraint
#' @param constraint_fn constraint function
#' @param constraint_grad_fn gradient of constraint fn
#' @param tolerance convergence tolerance. Once no element of B changes by 
#' this amount or more in an iteration, convergence is declared.
#' @param maxit maximum iterations
#' @param verbose shout at you?
#' @param trackB track value of beta across iterations and return?
#' @param starting_stepsize stepsize at which backtracking line search begins. Default is 0.5.
#' @param do_shift perform optimization of row-wise location of B in between rounds of coordinate descent?
#' Default is TRUE.
#' 
#' @returns A list containing elements `B`, `k_constr`, `j_constr`, `niter`
#' `converged`, `dll_dB`, and `Bs`. `B` is a matrix containing parameter estimates 
#' under the null (obtained by maximum likelihood on augmented observations Y),
#' 'k_constr', and 'j_constr' give row and column indexes of the parameter 
#' fixed to be equal to the constraint function g() under the null. `niter` is a 
#' scalar giving total number of outer iterations used to fit the null model.
#' 'converged' is an indicator that convergence was reached (TRUE) -- if optimization
#' stopped because the iteration limit was reached, `convergence` will have value FALSE. 
#' `dll_dB` is the derivative of the Poisson log likelihood with respect to the 
#' parameters in the reparametrized null model, and 
#' 'Bs' is a data frame containing values of B by iteration if trackB was set 
#' equal to TRUE (otherwise it contains a NULL value).
#' 
fit_null_repar <- function(B,
                           Y, 
                           X, 
                           k_constr, 
                           j_constr, 
                           j_ref, 
                           constraint_fn, 
                           constraint_grad_fn, 
                           tolerance = 1e-3,
                           maxit = 100, 
                           verbose = FALSE, 
                           trackB = FALSE,
                           starting_stepsize = 0.5,
                           method = "fisher",
                           do_shift = TRUE,
                           shift_tolerance = 1e-3
) {
  
  J <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  
  if(j_ref == j_constr){
    stop("Taxon j chosen for convenience constraint may not be null taxon.")
  }
  
  if(j_ref != J){
    stop("fit_null_repar has not yet been tested with j_ref != J.")
  }
  
  #change id. constr. by subtracting approp. col of B from others
  for(k in 1:p) {
    B[k,] <- B[k,] - B[k,j_ref]
  }
  
  # impose null constraint -- only valid for permutation-invariant constraint functions
  B[k_constr,j_constr] <- constraint_fn(B[k_constr,-j_constr])
  
  #update z
  z <- update_z(X = X, Y = Y, B = B)
  
  if(trackB){
    Bs <- vector(mode = "list",length = maxit + 1)
    Bs[[1]] <- B
  } else{
    Bs <- NULL
  }
  
  #evaluate log likelihood
  log_means <- X%*%B + matrix(z, ncol = 1)%*%matrix(1,ncol = J,nrow = 1)
  ll <- sum(Y*log_means - exp(log_means))
  lls <- ll
  
  prev_ll <- -Inf
  Bs <- list(B)
  
  #loop over j except j_constr and j_ref
  loop_j <- setdiff(1:J,c(j_constr,j_ref))
  
  #initiate prev_B at an arbitrary value that will not trigger convergence
  stop_iteration <- FALSE
  
  iter <- 1
  while(!stop_iteration){
    #store previous value of B
    prev_B <- B
    
    if(do_shift){
      lil_ll_rep <- function(shifts){
        temp_B <- B
        for(k in 1:p){
          temp_B[k,-j_ref] <-  temp_B[k,-j_ref] + shifts[k]
        }
        temp_B[k_constr,j_constr] <- constraint_fn( temp_B[k_constr,-j_constr])
        temp_z <- update_z(X = X, Y = Y, B = temp_B)
        log_means <- X%*%temp_B + matrix(temp_z, ncol = 1)%*%matrix(1,ncol = J,nrow = 1)
        
        return(-sum(Y*log_means - exp(log_means)))
      }
      
      shifty <- optim(rep(0,p),lil_ll_rep,
                      method = "BFGS")
     
      shift <- shifty$par
      
      if(max(abs(shift))<shift_tolerance ){
        do_shift <- FALSE
      }
    }
    
    if(do_shift){
      # print(signif(shift,2))
      
      for(k in 1:p){
        B[k,-j_ref] <-  B[k,-j_ref] + shift[k]
      }
      B[k_constr,j_constr] <- constraint_fn(B[k_constr,-j_constr])
      
      z <- update_z(X = X, Y = Y, B = B)
    }
    
    for(j in loop_j){
      # print(j)
      #compute update for B_j and unconstrained elements of B_j_constr
      update <- null_repar_one(B = B,
                               j = j,
                               z = z,
                               Y = Y,
                               X = X,
                               j_constr = j_constr,
                               k_constr = k_constr,
                               constraint_fn = constraint_fn,
                               constraint_grad_fn = constraint_grad_fn,
                               starting_stepsize = starting_stepsize,
                               method = method)
      
      #update B_j
      B[,j] <-  B[,j]  + update[1:p] 
      
      #update unconstrained elements of B_k_constr
      B[-k_constr,j_constr] <- B[-k_constr,j_constr] + update[(1:(p - 1))+p] 
      
      #update B_k_constr_j_constr to satisfy null
      B[k_constr,j_constr] <- constraint_fn(B[k_constr,-j_constr])
    }
    
    #update z
    z <- update_z(X = X, Y = Y, B = B)
    
    # compute ll
    log_means <- X%*%B + matrix(z, ncol = 1)%*%matrix(1,ncol = J,nrow = 1)
    prev_ll <- ll
    ll <- sum(Y*log_means - exp(log_means))
    if(verbose){
    message(round(ll,1))
    message(signif(ll - prev_ll,3))
    }
    
    lls <- c(lls,ll)
    Bs <- c(Bs,list(B))
    
    deriv_mat <- do.call(cbind,lapply(1:J,
                                      function(j) t(X)%*%(Y[,j,drop = FALSE] - exp(log_means[,j,drop = FALSE]))))
  
    deriv_mat[k_constr,-j_constr] <-  deriv_mat[k_constr,-j_constr]  + deriv_mat[k_constr,j_constr]*constraint_grad_fn(B[k_constr,-j_constr])
    deriv_mat[k_constr,j_constr] <- 0
    deriv_mat[,j_ref] <- 0
    
    if(verbose){
    message("Norm of derivative: ",round(sqrt(sum(deriv_mat^2)),2))
    if(iter>1){
      old_dl <- deriv_long
    }
    deriv_long <- do.call(c,lapply(1:p,function(k) deriv_mat[k,]))
    
    # 
    # if(iter ==1){
    #   plot(asinh(deriv_long),col = "red",pch = 20)
    # } else{
    #   points(asinh(old_dl), col = "grey",pch = 20)
    #   points(asinh(deriv_long),col = "red",pch = 20)
    # }
    }
    iter <- iter + 1
    
    if(iter > maxit){
      stop_iteration <- TRUE
      converged <- FALSE
    }
    
    if(verbose){
    message("Maximum absolute difference in B: ",signif(max(abs(B - prev_B)),3))
      }
    # hist(as.numeric(log(abs(B - prev_B))/log(10)),
    #      breaks = 100,xlab = "Log 10 absolute difference")
    if(max(abs(B - prev_B))<tolerance){
      stop_iteration <- TRUE
      converged <- TRUE
    }
    
    if(trackB){
      Bs[[iter]] <- B
    }
  }
  
  
  deriv_mat[k_constr,j_constr] <- NA
  
  return(list("B" = B,
              "k_constr" = k_constr,
              "j_constr" = j_constr,
              "niter" = iter - 1,
              "converged" = converged,
              "dll_dB" = deriv_mat,
              "Bs" = Bs))
}

