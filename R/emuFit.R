#' Fit radEmu model
#'
#' @param Y an n x J matrix of nonnegative observations
#' @param formula a one-sided formula specifying the form of the mean model to be fit
#' @param data an n x p data frame containing variables given in \code{formula}
#' @param B starting value of coefficient matrix (p x J). If not provided,
#' B will be initiated as a zero matrix.
#' @param test_kj a data frame whose rows give coordinates (in category j and
#' covariate k) of elements of B to construct hypothesis tests for. If \code{test_kj}
#' is not provided, all elements of B save the intercept row will be tested.
#' @param verbose provide updates as model is being fitted? Defaults to TRUE.
#' @param tolerance tolerance on squared norm of gradient of likelihood to be used
#' as a stopping condition in model fitting
#' @param maxit maximum number of coordinate descent cycles to perform before
#' exiting optimization
#' @return A p x J matrix containing regression coefficients (under constraint
#' g(B_k) = 0)
#' @param constraint_fn function g defining a constraint on rows of B; g(B_k) = 0
#' for rows k = 1, ..., p of B. Default function is a smoothed median (minimizer of
#' pseudohuber loss).
#' @param constraint_grad_fn derivative of constraint_fn with respect to its
#' arguments (i.e., elements of a row of B)
#' @param constraint_param If pseudohuber centering is used (this is the default),
#' parameter controlling relative weighting of elements closer and further from center.
#' (Limit as \code{constraint_param} approaches infinity is the mean; as this parameter approaches zero,
#' the minimizer of the pseudo-Huber loss approaches the median.)
#' 
#' @importFrom stats cov median model.matrix optim pchisq qnorm weighted.mean
#' @import Matrix
#' @import MASS
#' 
#' @export
#' 
emuFit <- function(Y,
                   formula,
                   data,
                   penalize = TRUE,
                   B = NULL,
                   test_kj = NULL,
                   verbose = TRUE,
                   tolerance = 1e-2,
                   maxit = 500,
                   tol_factor = 1e-2,
                   min_tol = 1e-2,
                   constraint_fn = pseudohuber_center,
                   constraint_grad_fn = dpseudohuber_center_dx,
                   constraint_param = 1,
                   alpha = 0.05,
                   return_wald_p = FALSE,
                   run_score_tests = TRUE,
                   rho_init = 1,
                   rho_scaling = 2,
                   gap_tolerance = 1e-6
){
  X <- model.matrix(formula,data)
  
  n <- nrow(Y)
  J <- ncol(Y)
  p <- ncol(X)
  
  
  if(all.equal(constraint_fn,pseudohuber_center)==TRUE){
    if(verbose){
      message("Centering rows of B with pseudo-Huber smoothed median with smoothing parameter ",constraint_param,".")
    }
  }
  
  
  if(constraint_param != 1){
    if(all.equal(constraint_fn,pseudohuber_center)==TRUE){
      constraint_fn <- (function(x) pseudohuber_center(x,d = constraint_param))
      constraint_grad_fn <- (function(x) dpseudohuber_center_dx(x,d = constraint_param))
      
    }
    else{
      warning("Argument constraint_param is currently only supported for centering with
pseudohuber_center() function; constraint_param input is otherwise ignored. Please directly
feed your choice of constraint function, including any necessary parameters, to constraint_fn argument
and the corresponding gradient function to constraint_grad_fn.")
    }
  }
  
  
  
  
  if(penalize){
    #fit penalized_model
    fitted_model <-
      emuFit_micro_penalized(X = X,
                             Y = Y,
                             B = B,
                             constraint_fn = constraint_fn,
                             maxit = maxit,
                             tolerance = tolerance,
                             verbose = verbose)
    Y_test <- fitted_model$Y_augmented
    fitted_B <- fitted_model$B
  } else{
    #fit ML model
    fitted_model <-
      emuFit_micro(X = X,
                   Y = Y,
                   B = B,
                   constraint_fn = NULL,
                   maxit = maxit,
                   tolerance = tolerance)
    fitted_B <- fitted_model
    Y_test <- Y
  }
  
  
  if(is.null(test_kj)){
    test_kj <- expand.grid(1:J, 2:p)
    test_kj <- data.frame(j = test_kj[,1],
                          k = test_kj[,2])
  }
  
  ntests <- nrow(test_kj)
  
  test_kj$estimate <- do.call(c, lapply(2:p, function(k) fitted_B[k,]))
  test_kj$lower <- NA
  test_kj$upper <- NA
  test_kj$pval <- test_kj$score_stat <- NA
  
  test_kj <- micro_wald(Y = Y,
                        X = X,
                        B = fitted_B,
                        test_kj = test_kj,
                        constraint_fn = constraint_fn,
                        constraint_grad_fn = constraint_grad_fn,
                        nominal_coverage = 1 - alpha)
  
  if(return_wald_p){
    test_kj$wald_p <- test_kj$p
  }
  test_kj <- test_kj[,colnames(test_kj) != "p"]
  
  
  
  for(test_ind in 1:ntests){
    if(run_score_tests){
      if(verbose){
        print(paste("Running score test ", test_ind, " of ", ntests," (row of B k = ", test_kj$k[test_ind],"; column of B j = ",
                    test_kj$j[test_ind],").",sep = ""))
      }
      
      H <- matrix(0,nrow = p, ncol = J)
      for(j in 1:J){
        H[null_k,j] <- constraint_grad_fn(fitted_B[test_kj$k[test_ind],])[j]
      }
      H[test_kj$k[test_ind],test_kj$j[test_ind]] <-  H[test_kj$k[test_ind],test_kj$j[test_ind]] - 1
      H_cup <- B_cup_from_B(H)
      H_cup[] <- as.numeric(H_cup[]!=0)
      
      
      
      
      
      test_result <- micro_score_test(Y = Y_test,
                                      X = X,
                                      B = fitted_B,
                                      constraint_fn = constraint_fn,
                                      constraint_grad_fn = constraint_grad_fn,
                                      null_k = test_kj$k[test_ind],
                                      null_j = test_kj$j[test_ind],
                                      rho_init = rho_init,
                                      rho_scaling = rho_scaling,
                                      gap_tolerance = gap_tolerance,
                                      maxit = 500,
                                      tolerance = tolerance,
                                      # step_ratio = 0.1,
                                      verbose = verbose)
    }
    #
    #       message("David to test: pseudohuber derivatives with at steps where ll jumps around;
    # augmented ll derivatives at similar steps; info positive definite (which it should always be); probably more.")
    
    test_kj$estimate[test_ind] <- fitted_B[test_kj$k[test_ind],test_kj$j[test_ind]]
    
    if(run_score_tests){
      test_kj$pval[test_ind] <- test_result$pval
      test_kj$score_stat[test_ind] <- test_result$score_stat
    } else{
      test_kj$pval[test_ind] <- NA
      test_kj$score_stat[test_ind] <- NA
    }
    
    
  }
  
  k_to_covariates <- data.frame(k = 1:p,
                                covariate = colnames(X))
  
  test_kj$covariate <- do.call(c,
                               lapply(test_kj$k,
                                      function(d) k_to_covariates$covariate[
                                        k_to_covariates$k ==d
                                      ]))
  
  if(!is.null(colnames(Y))){
    j_to_categories <- data.frame(j = 1:J,
                                  category = colnames(Y))
    test_kj$category <- do.call(c,
                                lapply(test_kj$j,
                                       function(d) j_to_categories$category[
                                         j_to_categories$j ==d
                                       ]))
  } else{
    test_kj$category <- NA
  }
  
  test_kj <-
    cbind(data.frame(covariate = test_kj$covariate,
                     category = test_kj$category,
                     category_num = test_kj$j),
          test_kj[,c("estimate","se","lower","upper","score_stat","pval")])
  
  
  # message("David to add: rho_init, rho_scaling, gap_tolerance")
  # message("David also to add: labeling of output by input (taxa) and labeling of input by output (covariates)")
  # message("David as well adding: option to not perform score tests")
  
  return(test_kj)
}


