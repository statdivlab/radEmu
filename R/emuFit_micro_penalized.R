
#' Fit radEmu model with Firth penalty
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
#' @param max_step numeric: maximum sup-norm for proposed update steps
#' @param verbose logical: report information about progress of optimization? Default is TRUE.
#' @return A p x J matrix containing regression coefficients (under constraint
#' g(B_k) = 0)
#' 
emuFit_micro_penalized <-
  function(X,
           Y,
           B = NULL,
           constraint_fn = NULL,
           maxit = 500,
           tolerance = 1e-5,
           max_step = 0.5,
           verbose = TRUE
           ){

J <- ncol(Y)
p <- ncol(X)
n <- nrow(Y)
X_tilde <- X_cup_from_X(X,J)
Y_augmented <- Y
fitted_model <- NULL
converged <- FALSE
counter <- 0
while(!converged){
  
  if(counter ==0){
    Y_augmented <- Y + 1
  }

fitted_model <- emuFit_micro(X,
                            Y_augmented,
                            B = fitted_model,
                            constraint_fn = constraint_fn,
                            maxit = maxit,
                            max_stepsize = max_step,
                            tolerance = tolerance,
                            verbose = verbose)

Y_augmented <-
  update_data(Y = Y,
              X_tilde = X_tilde,
              B = fitted_model,
              p = p,
              n = n,
              J = J)

deriv <- dpll_dB_cup(X,Y_augmented,fitted_model)
deriv_norm  = sqrt(sum(deriv^2))

if(deriv_norm<tolerance){
  converged <- TRUE
}

if(counter>maxit){
  converged <- TRUE
}
counter <- counter + 1
}

return(list("Y_augmented" = Y_augmented,
            "B" = fitted_model))

  }
