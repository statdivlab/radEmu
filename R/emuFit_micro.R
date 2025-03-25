
#' Fit radEmu model. Called by emuFit
#'
#' @param X a p x J design matrix
#' @param Y an n x p matrix of nonnegative observations
#' @param B starting value of coefficient matrix (p x J)
#' @param constraint_fn function g defining constraint on rows of B; g(B_k) = 0
#' for rows k = 1, ..., p of B.
#' @param maxit maximum number of coordinate descent cycles to perform before
#' exiting optimization
#' @param tolerance tolerance on improvement in log likelihood at which to
#' exit optimization
#' @param verbose logical: print information about optimization progress? Default is TRUE.
#' @param warm_start logical: begin from "warm start" obtained from linear regression on
#' transformed counts? Default is TRUE.
#' @param c1 numeric: value of constant in Armijo condition checked in backtracking line search
#' @param max_stepsize numeric: maximum sup-norm value of proposed step. Default is 0.5.
#' @param max_abs_B numeric: maximum value elements of B are allowed to take in absolute value.
#' Helps prevent optimization failure in larger problems. Defaults to 50.
#' @param use_working_constraint logical: set a column of B equal to
#' zero within optimization. Default is TRUE.
#' @param j_ref If use_working_constraint is TRUE, column index of column of B
#' to set to zero. Default is NULL, in which case this column is chosen to
#' maximize the number of nonzero entries of Y_j_ref.
#' @param optimize_rows If use_working_constraint is TRUE, update overall location of 
#' rows of B relative to column constrained to equal zero under working constraint before
#' iterating through updates to columns of B individually. Default is TRUE.
#' @return A p x J matrix containing regression coefficients (under constraint
#' g(B_k) = 0)
#'
#' @export
emuFit_micro <-
  function(X,
           Y,
           B = NULL,
           constraint_fn = NULL,
           maxit = 250,
           tolerance = 1e-5,
           verbose = TRUE,
           warm_start = TRUE,
           c1 = 1e-4,
           max_stepsize = 0.5,
           max_abs_B = 50,
           use_working_constraint = TRUE,
           j_ref = NULL,
           optimize_rows = TRUE){
    #extract dimensions
    n <- nrow(Y)
    J <- ncol(Y)
    p <- ncol(X)


    # assign constraint if necessary
    if(is.null(constraint_fn)){
      constraint_fn <- rep(list(function(x) pseudohuber_center(x,0.1)), p)
      constraint_grad_fn <- rep(list(function(x) dpseudohuber_center_dx(x,0.1)), p)
    }


    #populate B if not provided
    if(is.null(B)){
      if(!warm_start){
      B <- matrix(0,
                  nrow = p,
                 ncol = J)
      } else{

        Y_start <- Y + 0.01*mean(Y)
        for(i in 1:nrow(Y_start)){
          Y_start[i,] <- log(Y_start[i,])
          Y_start[i,] <- Y_start[i,] - mean(Y_start[i,])
        }
        B <- matrix(nrow = p,ncol = J)
        
        ##Checking if design matrix is rank-deficient and soln_mat can be found
        if (qr(t(X)%*%X)$rank < ncol(X)){
          stop("Design matrix X inputted for the model is rank-deficient, preventing proper model fitting.
  This might be due to multicollinearity, overparameterization, or redundant factor levels included in covariates.
  Consider removing highly correlated covariates or adjusting factor levels to ensure a full-rank design. \n")
        }
        
        soln_mat <- qr.solve(t(X)%*%X,t(X))
        for(j in 1:J){
          B[,j] <- as.numeric(as.matrix(soln_mat%*%Y_start[,j,drop = FALSE]))
        }
      }
    }



    loop_js <- 1:J
    if(use_working_constraint){
      if(is.null(j_ref)){
      detection <- colSums(Y>0)

      which_to_omit <- which(detection == max(detection))
      if(length(which_to_omit)>1){
        totals <- colSums(Y[,which_to_omit])
        which_to_omit <- which_to_omit[which.max(totals)]
      }
      } else{
        which_to_omit <- j_ref
      }

      loop_js <- loop_js[-which_to_omit]

      for(k in 1:p){
        B[k,] <- B[k,] - B[k,which_to_omit]
      }
    } else{
      for(k in 1:p){
        B[k,] <- B[k,] - constraint_fn[[k]](B[k,])
      }
    }




    #create object to store log likelihoods at each iteration in
    lls <- numeric(maxit)

    #initiate z
    z <- update_z(Y,X,B)

    #compute initial mean
    log_mean <- X%*%B +
      matrix(z,ncol = 1)%*%matrix(1,ncol = J, nrow = 1)

    lls[1] <- sum(Y*log_mean - exp(log_mean))

    converged <- FALSE

    if(optimize_rows&use_working_constraint){
      
      eps_outcome <- rowSums(Y[,-which_to_omit,drop = FALSE])
      
    }


    iter <- 1
    deriv_norm <- Inf
    B_diff <- Inf
    iter <- 1
    while(!converged){
      old_B <- B
      if(optimize_rows &use_working_constraint){
        
        epsilon_step <- rep(10,p)
        while(max(abs(epsilon_step))>0.01){

        log_mean <- X%*%B +
          matrix(z,ncol = 1)%*%matrix(1,ncol = J, nrow = 1)
        eps_offset =  log(rowSums(exp(log_mean[,-which_to_omit,drop= FALSE])))
        epsilon_step <- suppressWarnings(stats::glm(eps_outcome~X-1,family = "poisson",offset =eps_offset,
                       control = list(maxit = 3))$coef)

            for(k in 1:p){
          B[k,-which_to_omit] <-  B[k,-which_to_omit] + epsilon_step[k]
            }
          z <- update_z(Y,X,B)
        }}
        

      for(j in loop_js){
      
        update <- micro_fisher(X,Yj = Y[,j,drop= FALSE],Bj = B[,j,drop = FALSE],z,
                               stepsize = max_stepsize,
                               c1 = c1)
        
        B[,j] <- B[,j] + update

        if(!use_working_constraint){
        for(k in 1:p){
          B[k,] <- B[k,] - constraint_fn[[k]](B[k,])
        }}



      }

      # keep any elements of B from flying off to infinity
      B[B < -max_abs_B] <- -max_abs_B
      B[B > max_abs_B] <- max_abs_B
      z <- update_z(Y = Y,X = X,B = B)



      log_mean <- X%*%B +
        matrix(z,ncol = 1)%*%matrix(1,ncol = J, nrow = 1)
      lls[iter  +1] <- sum(Y*log_mean - exp(log_mean))

        deriv <- do.call(c,lapply(1:J,
                                  function(j) crossprod(X,Y[,j,drop = FALSE] - exp(log_mean[,j,drop = FALSE]))))

        deriv_norm  = sqrt(sum(deriv^2))

        B_diff <- max(abs(B - old_B)[abs(B)<(0.5*max_abs_B)])

      if(verbose){
      message("Scaled norm of derivative ",signif(sqrt((1/(n*J))*deriv_norm),4) )
        if(iter>1){
          message("Max absolute difference in elements of B since last iteration: ",signif(B_diff,3))
        }
        }

      # if(deriv_norm < tolerance){
      #   converged <- TRUE
      # }
        if(B_diff < tolerance){
          converged <- TRUE
        }


      if(iter > maxit){
        converged <- TRUE
        if(verbose){
        message("Iteration limit reached; exiting optimization.")
          }
      }
      iter <- iter + 1
    }

    if(use_working_constraint){
      for(k in 1:p){
        B[k,] <- B[k,] - constraint_fn[[k]](B[k,])
      }}

    z <- update_z(Y = Y,X = X,B = B)

    log_mean <- X%*%B +
      matrix(z,ncol = 1)%*%matrix(1,ncol = J, nrow = 1)
    lls[iter  +1] <- sum(Y*log_mean - exp(log_mean))

    deriv <- do.call(c,lapply(1:J,
                              function(j) crossprod(X,Y[,j,drop = FALSE] - exp(log_mean[,j,drop = FALSE]))))

    deriv_norm  = sqrt(sum(deriv^2))
    # print(deriv_norm)

      return(B)
  }

