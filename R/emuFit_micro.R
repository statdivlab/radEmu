
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
           tolerance = 1e-4,
           verbose = TRUE,
           warm_start = TRUE,
           c1 = 1e-4,
           max_stepsize = 0.5){
    #extract dimensions
    n <- nrow(Y)
    J <- ncol(Y)
    p <- ncol(X)


    # assign constraint if necessary
    if(is.null(constraint_fn)){
      constraint_fn <- pseudohuber_center
      constraint_grad_fn <- dpseudohuber_center_dx
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
        soln_mat <- qr.solve(t(X)%*%X,t(X))
        for(j in 1:J){
          B[,j] <- soln_mat%*%Y_start[,j,drop = FALSE]
        }
      }
    }

    for(k in 1:p){
      B[k,] <- B[k,] - constraint_fn(B[k,])
    }




    #create object to store log likelihoods at each iteration in
    lls <- numeric(maxit)

    #initiate z
    z <- update_z_no_wts(Y,X,B)

    #compute initial mean
    log_mean <- X%*%B +
      matrix(z,ncol = 1)%*%matrix(1,ncol = J, nrow = 1)

    lls[1] <- sum(Y*log_mean - exp(log_mean))

    converged <- FALSE


    iter <- 1
    deriv_norm <- Inf
    while(!converged){
      for(j in 1:J){
        # print(j)
        llj <- function(x){
          Bj <- matrix(x,ncol = 1)
          log_mean <- X%*%Bj + z
          return(sum(Y[,j]*log_mean - exp(log_mean)))
        }
        grj <- function(x){
          Bj <- matrix(x,ncol = 1)
          log_mean <- X%*%Bj + z
          return(-1*apply(diag(as.numeric((Y[,j] - exp(log_mean))))%*%X,
                       2,sum))
        }

        B[,j] <- B[,j] + micro_fisher(X,Yj = Y[,j,drop= FALSE],Bj = B[,j,drop = FALSE],z,
                                      stepsize = max_stepsize,
                                      c1 = c1)


        for(k in 1:p){
          B[k,] <- B[k,] - constraint_fn(B[k,])
        }

        z <- update_z_no_wts(Y = Y,X = X,B = B)

      }

      log_mean <- X%*%B +
        matrix(z,ncol = 1)%*%matrix(1,ncol = J, nrow = 1)
      lls[iter  +1] <- sum(Y*log_mean - exp(log_mean))

      if(J>100){
        deriv <- do.call(c,lapply(1:J,
                                  function(j) crossprod(X,Y[,j,drop = FALSE] - exp(log_mean[,j,drop = FALSE]))))
        deriv_norm  = sqrt(sum(deriv^2))
      } else{
      deriv <- dpll_dB_cup(X,Y,B)
      deriv_norm  = sqrt(sum(deriv^2))
      }

      if(verbose){
      message("Norm of derivative ",signif(deriv_norm,4) )
        }

      if(deriv_norm < tolerance){
        converged <- TRUE
      }


      if(iter > maxit){
        converged <- TRUE
        message("Iteration limit reached; exiting optimization.")
      }
      iter <- iter + 1
    }

      return(B)
  }
