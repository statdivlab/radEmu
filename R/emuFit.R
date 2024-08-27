#' Fit radEmu model
#'
#' @param Y an n x J matrix or dataframe of nonnegative observations, or a phyloseq object containing an otu table and sample data.
#' @param X an n x p matrix or dataframe of covariates (optional)
#' @param formula a one-sided formula specifying the form of the mean model to be fit
#' @param data an n x p data frame containing variables given in \code{formula}
#' @param cluster a vector giving cluster membership for each row of Y to be used in computing 
#' GEE test statistics. Default is NULL, in which case rows of Y are treated as independent.
#' @param penalize logical: should Firth penalty be used in fitting model? Default is TRUE.
#' @param B starting value of coefficient matrix (p x J). If not provided,
#' B will be initiated as a zero matrix.
#' @param B_null_list list of starting values of coefficient matrix (p x J) for null estimation. This should either 
#' be a list with the same length as \code{test_kj}. If you only want to provide starting values for some tests,
#' include the other elements of the list as \code{NULL}.
#' @param fitted_model a fitted model produced by a separate call to emuFit; to
#' be provided if score tests are to be run without refitting the full unrestricted model.
#' Default is NULL.
#' @param refit logical: if B or fitted_model is provided, should full model be fit (TRUE) or
#' should fitting step be skipped (FALSE), e.g., if score tests are to be run on an already
#' fitted model. Default is TRUE.
#' @param test_kj a data frame whose rows give coordinates (in category j and
#' covariate k) of elements of B to construct hypothesis tests for. If \code{test_kj}
#' is not provided, all elements of B save the intercept row will be tested.
#' @param alpha nominal type 1 error level to be used to construct confidence intervals. Default is 0.05
#' (corresponding to 95% confidence intervals)
#' @param return_wald_p logical: return p-values from Wald tests? Default is FALSE.
#' @param compute_cis logical: compute and return Wald CIs? Default is TRUE.
#' @param run_score_tests logical: perform robust score testing? Default is TRUE.
#' @param use_fullmodel_info logical: TODO? Default is FALSE.
#' @param use_fullmodel_cov logical: use information matrix and empirical score covariance
#' computed for full model fit? Defaults to FALSE, in which case these quantities are
#' recomputed for each null model fit for score testing.
#' @param use_both_cov logical: should score tests be run using information and
#' empirical score covariance evaluated both under the null and full models?
#' Used in simulations
#' @param constraint_fn function g defining a constraint on rows of B; g(B_k) = 0
#' for rows k = 1, ..., p of B. Default function is a smoothed median (minimizer of
#' pseudohuber loss). If a number is provided a single category constraint will be used
#' with the provided category as a reference category. 
#' @param constraint_grad_fn derivative of constraint_fn with respect to its
#' arguments (i.e., elements of a row of B)
#' @param constraint_param If pseudohuber centering is used (this is the default),
#' parameter controlling relative weighting of elements closer and further from center.
#' (Limit as \code{constraint_param} approaches infinity is the mean; as this parameter approaches zero,
#' the minimizer of the pseudo-Huber loss approaches the median.)
#' @param verbose provide updates as model is being fitted? Defaults to FALSE.
#' @param tolerance tolerance for stopping criterion in full model fitting; once
#' no element of B is updated by more than this value in a single step, we exit
#' optimization. Defaults to 1e-3.
#' @param rho_init numeric: value at which to initiate rho parameter in augmented Lagrangian
#' algorithm. Default is 1.
#' @param tau numeric: value to scale rho by in each iteration of augmented Lagrangian
#' algorithm that does not move estimate toward zero sufficiently. Default is 2.
#' @param kappa numeric: value between 0 and 1 that determines the cutoff on the ratio
#' of current distance from feasibility over distance in last iteration triggering
#' scaling of rho. If this ratio is above kappa, rho is scaled by tau to encourage
#' estimate to move toward feasibility.
#' @param constraint_tol numeric: constraint tolerance for fits under null hypotheses
#' (tested element of B must be equal to constraint function to within this tolerance for
#' a fit to be accepted as a solution to constrained optimization problem). Default is 1e-5.
#' @param maxit maximum number of outer iterations of augmented lagrangian algorithm to perform before
#' exiting optimization. Default is 1000.
#' @param inner_maxit maximum number of coordinate descent passes through columns of B to make within each
#' outer iteration of augmented lagrangian algorithm before exiting inner loop
#' @param max_step maximum stepsize; update directions computed during optimization
#' will be rescaled if a step in any parameter exceeds this value. Defaults to 0.5.
#' @param B_null_tol numeric: convergence tolerance for null model fits for score testing (if max of absolute difference in
#' B across outer iterations is below this threshold, we declare convergence).
#' Default is 0.01.
#' @param inner_tol numeric: convergence tolerance for augmented Lagrangian subproblems
#' within null model fitting. Default is 1.
#' @param ntries numeric: how many times should optimization be tried in null
#' models where at least one optimization attempt fails? Default is 4.
#' @param c1 numeric: parameter for Armijo line search. Default is 1e-4.
#' @param trackB logical: should values of B be recorded across optimization
#' iterations and be returned? Primarily used for debugging. Default is FALSE.
#' @param return_nullB logical: should values of B under null hypothesis be returned. Primarily used for debugging. Default is FALSE. 
#' @param return_both_score_pvals logical: should score p-values be returned using both
#' information matrix computed from full model fit and from null model fits? Default is
#' FALSE. This parameter is used for simulations - in any applied analysis, type of
#' p-value to be used should be chosen before conducting tests.
#' @param remove_zero_comparison_pvals Should score p-values be replaced with NA for zero-comparison parameters? These parameters occur 
#' for categorical covariates with three or more levels, and represent parameters that compare a covariate level to the reference level for
#' a category in which the comparison level and reference level both have 0 counts in all samples. These parameters can have misleadingly 
#' small p-values and are not thought to have scientifically interesting signals. We recommend removing them before analyzing data further. 
#' If TRUE, all zero-comparison parameter p-values will be set to NA. If FALSE no zero-comparison parameter p-values will be set to NA.
#' If a value between 0 and 1, all zero-comparison p-values below the value will be set to NA. 
#' Default is \code{0.01}. 
#' 
#' @return A list containing elements 'coef', 'B', 'penalized', 'Y_augmented',
#' 'z_hat', 'I', 'Dy', and 'score_test_hyperparams' if score tests are run.  
#' Parameter estimates by covariate and outcome category (e.g., taxon for microbiome data),
#' as well as optionally confidence intervals and p-values, are contained in 'coef'.
#' Any robust score statistics and score test p-values are also included in 'coef'. 
#' If there are any zero-comparison parameters in the model, a column 'zero_comparison'
#' is also included, which is TRUE for any parameters that compare the level of a categorical
#' covariate to a reference level for a category with only zero counts for both the comparison
#' level and the reference level. This check is currently implemented for an arbitrary design matrix
#' generated using the \code{formula} and \code{data} arguments, and for a design matrix with no more
#' than one categorical covariate if the design matrix \code{X} is input directly.
#' 'B' contains parameter estimates in matrix format (rows indexing covariates and
#' columns indexing outcome category / taxon). 
#' 'penalized' is equal to TRUE f Firth penalty is used in estimation (default) and FALSE otherwise. 
#' 'z_hat' returns the nuisance parameters calculated in Equation 7 of the radEmu manuscript,
#' corresponding to either 'Y_augmented' or 'Y' if the 'penalized' is equal to TRUE
#' or FALSE, respectively. 
#' I' and 'Dy' contain an information matrix and empirical score covariance matrix computed under the 
#' full model. 
#' 'score_test_hyperparams' contains parameters and hyperparameters related to estimation under the null,
#' including whether or not the algorithm converged, which can be helpful for debugging. 
#'
#' @importFrom stats cov median model.matrix optim pchisq qnorm weighted.mean
#' @import Matrix
#' @import MASS
#'
#' @export
#'
emuFit <- function(Y,
                   X = NULL,
                   formula = NULL,
                   data = NULL,
                   cluster = NULL,
                   penalize = TRUE,
                   B = NULL,
                   B_null_list = NULL,
                   fitted_model = NULL,
                   refit = TRUE,
                   test_kj = NULL,
                   alpha = 0.05,
                   return_wald_p = FALSE,
                   compute_cis = TRUE,
                   run_score_tests = TRUE,
                   use_fullmodel_info = FALSE,
                   use_fullmodel_cov = FALSE,
                   use_both_cov = FALSE,
                   constraint_fn = pseudohuber_center,
                   constraint_grad_fn = dpseudohuber_center_dx,
                   constraint_param = 0.1,
                   verbose = FALSE,
                   tolerance = 1e-4,
                   B_null_tol = 1e-3,
                   rho_init = 1,
                   inner_tol = 1,
                   ntries = 4,
                   tau = 2,
                   kappa = 0.8,
                   constraint_tol = 1e-5,
                   c1 = 1e-4,
                   maxit = 1000,
                   inner_maxit = 25,
                   max_step = 1,
                   trackB = FALSE,
                   return_nullB = FALSE,
                   return_both_score_pvals = FALSE,
                   remove_zero_comparison_pvals = 0.01) {
  
  # Record call
  call <- match.call(expand.dots = FALSE)
  
  # check if Y is a phyloseq object
  if ("phyloseq" %in% class(Y)) {
    if (requireNamespace("phyloseq", quietly = TRUE)) {
      if (is.null(formula)) {
        stop("If Y is a `phyloseq` object, make sure to include the formula argument.")
      } else {
        data <- data.frame(phyloseq::sample_data(Y))
        X <- model.matrix(formula, data)
        taxa_are_rows <- Y@otu_table@taxa_are_rows
        Y <- as.matrix(phyloseq::otu_table(Y))
        if (taxa_are_rows) {
          Y <- t(Y)
        }
      }
    } else {
      stop("You are trying to use a `phyloseq` data object or `phyloseq` helper function without having the `phyloseq` package installed. Please either install the package or use a standard data frame.")
    }
  } else if ("data.frame" %in% class(Y)) {
    Y <- as.matrix(Y)
    if (!is.numeric(Y)) {
      stop("Y is a data frame that cannot be coerced to a numeric matrix. Please fix and try again.")
    }
  }
  
  if (is.null(X)) {
    if (is.null(formula) | is.null(data)) {
      stop("If design matrix X not provided, both formula and data containing
covariates in formula must be provided.")
    }
    X <- model.matrix(formula, data)
  }
  if ("data.frame" %in% class(X)) {
    X <- as.matrix(X)
    if (!is.numeric(X)) {
      stop("X is a data frame that cannot be coerced to a numeric matrix. Please fix and try again.")
    }
  }
  
  # check that if X and Y have rownames, they match 
  if (!is.null(rownames(Y)) & !is.null(rownames(X))) {
    if (all.equal(rownames(Y), rownames(X)) != TRUE) {
      message("There is a different row ordering between covariate data and response data. Covariate data will be reordered to match response data.")
      X <- X[rownames(Y), ]
    }
  }
  
  if (min(rowSums(Y))==0) {
    stop("Some rows of Y consist entirely of zeroes, meaning that some samples
have no observations. These samples must be excluded before fitting model.")
  }
  
  if (min(colSums(Y)) == 0) {
    stop("Some columns of Y consist entirely of zeroes, meaning that some categories have zero counts for all samples. These
         categories must be excluded before fitting the model.")
  }
  
  #check that cluster is correctly type if provided
  if(!is.null(cluster)){
    if(length(cluster)!=nrow(Y)){
        stop("If provided as a vector, argument 'cluster' must have 
length equal to n (the number of rows in Y).")
      }
    if(length(unique(cluster)) == nrow(Y)){
      warning("Number of unique values in 'cluster' equal to number of rows of Y; 
ignoring argument 'cluster'.")
      cluster <- NULL
    }
  }
  
  # check that B_null_list is the correct length if provided 
  if (!is.null(B_null_list)) {
    if (length(B_null_list) != nrow(test_kj)) {
      warning("Length of 'B_null_list' is different than the number of tests specified in 'test_kj'. Ignoring object 'B_null_list'.")
      B_null_list <- NULL
    }
  }
  
  # check for valid argument remove_zero_comparison_pvals 
  if (remove_zero_comparison_pvals != TRUE & remove_zero_comparison_pvals != FALSE) {
    if (!(is.numeric(remove_zero_comparison_pvals) & remove_zero_comparison_pvals <= 1 &
          remove_zero_comparison_pvals >= 0)) {
      stop("Please set `remove_zero_comparison_pvals` to either TRUE, FALSE, or a numeric value between 0 and 1.")
    }
  }
  
  if (length(constraint_fn) == 1 & is.numeric(constraint_fn)) {
    constraint_cat <- constraint_fn
    constraint_fn <- function(x) {x[constraint_cat]}
    constraint_grad_fn <- function(x) {
      grad <- rep(0, length(x))
      grad[constraint_cat] <- 1
      return(grad)
    }
    constraint_param <- NA
  }
  
  n <- nrow(Y)
  J <- ncol(Y)
  p <- ncol(X)
  
  if (is.null(colnames(X))) {
    if (p > 1) {
      colnames(X) <- c("Intercept", paste0("covariate_", 1:(ncol(X) - 1)))
    } else {
      colnames(X) <- "Intercept"
    }
  }
  if (is.null(colnames(Y))) {
    colnames(Y) <- paste0("category_", 1:ncol(Y))
  }
  
  # check for zero-comparison parameters
  zero_comparison_res <- zero_comparison_check(X = X, Y = Y)
  
  X_cup <- X_cup_from_X(X,J)
  
  
  if (is.logical(all.equal(constraint_fn, pseudohuber_center))) {
    if (all.equal(constraint_fn, pseudohuber_center)) {
      if (verbose) message("Centering rows of B with pseudo-Huber smoothed median with smoothing parameter ", constraint_param, ".")
      
      stopifnot(!is.na(constraint_param))
      
      constraint_fn <- (function(x) pseudohuber_center(x, d = constraint_param))
      constraint_grad_fn <- (function(x) dpseudohuber_center_dx(x, d = constraint_param))
      
    } 
  } else {
    
    if (!is.na(constraint_param)) {
      
      warning("Argument constraint_param is currently only supported for centering with
pseudohuber_center() function; constraint_param input is otherwise ignored. Please directly
feed your choice of constraint function, including any necessary parameters, to constraint_fn argument
and the corresponding gradient function to constraint_grad_fn.")
    }
  } 
  
  #choose ref taxon for fitting constrained models / performing wald and score tests
  j_ref <- get_j_ref(Y)
  
  ########################################
  ########## fit model if needed
  ########################################
  
  if (refit) {

    
    if (penalize) {
      
      fitted_model <-
        emuFit_micro_penalized(X = X,
                               Y = Y,
                               B = B,
                               X_cup = X_cup,
                               constraint_fn = constraint_fn,
                               maxit = maxit,
                               max_step = max_step,
                               tolerance = tolerance,
                               verbose = verbose,
                               j_ref = j_ref)
      Y_test <- fitted_model$Y_augmented
      fitted_B <- fitted_model$B
      converged_estimates <- fitted_model$convergence
      if (!converged_estimates) {
        warning("Optimization to estimate parameters did not converge. Try running again with a larger value of 'maxit'.")
      }
    } else {
      
      fitted_model <-
        emuFit_micro(X = X,
                     Y = Y,
                     B = B,
                     constraint_fn = NULL,
                     maxit = maxit,
                     max_stepsize = max_step,
                     tolerance = tolerance,
                     j_ref = j_ref,
                     verbose = verbose)
      fitted_B <- fitted_model
      Y_test <- Y
    }
    
  } else {
    if (is.null(B) & is.null(fitted_model)) {
      stop("If refit is FALSE, B or fitted_model must be provided")
    }
    if (!is.null(fitted_model)) {
      fitted_B <- fitted_model$B
      if (penalize != fitted_model$penalized) {
        stop("Your argument to `penalize` does not match the `penalize` argument within your `fitted_model` object. Please use the `penalize` argument that matches the `penalized` return object within your `fitted_model`.")
      }
      if (penalize) {
        Y_test <- Y_augmented <- fitted_model$Y_augmented
      } else {
        Y_augmented <- NULL
        Y_test <- Y
      }
      if (!is.null(B)) {
        warning("B and fitted_model provided to emuFit; B ignored in favor of fitted_model.")
      }
    } else {
      fitted_B <- B
      if (penalize) {
        X_cup <- X_cup_from_X(X, J)
        G <- get_G_for_augmentations(X, J, n, X_cup)
        Y_test <- Y_augmented <- Y + 
          get_augmentations(X = X, G = G, Y = Y, B = fitted_B)
      } else {
        Y_augmented <- NULL
        Y_test <- Y
      }
    }
  }
  
  if (p == 1) {
    message("You are running an intercept-only model. In this model the intercept beta_0^j is an unidentifiable combination of intercept for category j and the detection efficiency of category j. Therefore this parameter is not interpretable.")
    
    coefficients <- expand.grid(1:J, 1)
    coefficients <- data.frame(j = coefficients[ , 1],
                               k = coefficients[ ,2])
    coefficients$estimate <- fitted_B[1, ]
    coefficients$lower <- NA
    coefficients$upper <- NA
    coefficients$pval <- coefficients$score_stat <- NA
  } else {
    coefficients <- expand.grid(1:J, 2:p)
    coefficients <- data.frame(j = coefficients[ , 1],
                               k = coefficients[ ,2])
    coefficients$estimate <- do.call(c, lapply(2:p, function(k) fitted_B[k,]))
    coefficients$lower <- NA
    coefficients$upper <- NA
    coefficients$pval <- coefficients$score_stat <- NA
  }
  
  if (compute_cis) {
    if (verbose) {
      message("Performing Wald tests and constructing CIs.")
    }
    
    just_wald_things <-  micro_wald(Y = Y_test,
                                    X = X,
                                    X_cup = X_cup,
                                    B = fitted_B,
                                    test_kj = coefficients,
                                    constraint_fn = constraint_fn,
                                    constraint_grad_fn = constraint_grad_fn,
                                    nominal_coverage = 1 - alpha,
                                    verbose = verbose,
                                    j_ref = j_ref,
                                    cluster = cluster)
    
    coefficients <- just_wald_things$coefficients
    if (use_fullmodel_cov) {
      I <- just_wald_things$I
      Dy <- just_wald_things$Dy
    } else {
      I <- NULL
      Dy <- NULL
    }
  } else {
    just_wald_things <- NULL
    I <- NULL
    Dy <- NULL
  }
  
  if (return_wald_p) {
    if (!compute_cis) {
      warning("Wald p-values cannot be returned if compute_cis = FALSE.")
    } else {
      coefficients$wald_p <- coefficients$pval}
  }
  
  coefficients$pval <- NA
  
  if (is.null(test_kj)) {
    test_kj <- coefficients
  }
  
  #add column for score p value from score stat using full model I, Dy if needed
  if (use_both_cov) {
    coefficients$score_fullcov_p <- NA
  }
  
  if (use_fullmodel_info) {
    indexes_to_remove = (j_ref - 1)*p + 1:p
    I_inv <- Matrix::solve(just_wald_things$I[-indexes_to_remove,-indexes_to_remove],  method = "cholmod_solve")
  } else {
    if(use_both_cov){
      I_inv <- Matrix::solve(just_wald_things$I[-indexes_to_remove,-indexes_to_remove],  method = "cholmod_solve")
    } else{
    I_inv <- NULL
    }
  }
  
  if (run_score_tests) {
    
    score_test_hyperparams <- data.frame(u = rep(NA, nrow(test_kj)),
                                         rho = NA,
                                         tau = NA,
                                         inner_maxit = NA,
                                         gap = NA,
                                         converged = NA)
    
    if (return_both_score_pvals) {
      colnames(coefficients)[colnames(coefficients) == "pval"] <-
        "score_pval_full_info"
      colnames(coefficients)[colnames(coefficients) == "score_stat"] <-
        "score_stat_full_info"
      coefficients$score_pval_null_info <- NA
      coefficients$score_stat_null_info <- NA
      
      if (!use_fullmodel_info) {
        stop("If return_both_score_pvals = TRUE, use_fullmodel_info must be TRUE as well.")
      }
    }
    
    if (return_nullB) {
      nullB_list <- vector(mode = "list", length = nrow(test_kj))
    }
    for(test_ind in 1:nrow(test_kj)) {
      
      if (verbose) {
        print(paste("Running score test ", test_ind, " of ", nrow(test_kj)," (row of B k = ", test_kj$k[test_ind], "; column of B j = ",
                    test_kj$j[test_ind],").",sep = ""))
      }
      
      B_to_use <- fitted_B
      if (!is.null(B_null_list)) {
        if (!is.null(B_null_list[[test_ind]])) {
          B_to_use <- B_null_list[[test_ind]]
          if (!(nrow(B_to_use) == nrow(fitted_B) & ncol(B_to_use) == ncol(fitted_B))) {
            warning("'B_null_list' contains objects that are not the correct dimension for 'B'. The 'B_null_list' argument will be ignored.")
            B_to_use <- fitted_B
          }
        }
      }
      
      test_result <- score_test(B = B_to_use, #B (MPLE or starting value if provided)
                                Y = Y_test, #Y (with augmentations)
                                X = X, #design matrix
                                X_cup = X_cup,
                                k_constr = test_kj$k[test_ind], #row index of B to constrain
                                j_constr = test_kj$j[test_ind],
                                constraint_fn = constraint_fn, #constraint function
                                constraint_grad_fn = constraint_grad_fn, #gradient of constraint fn
                                rho_init = rho_init,
                                tau = tau,
                                kappa = kappa,
                                B_tol = B_null_tol,
                                inner_tol = inner_tol,
                                constraint_tol = constraint_tol,
                                j_ref = j_ref,
                                c1 = c1,
                                maxit = maxit,
                                inner_maxit = inner_maxit,
                                ntries = ntries,
                                verbose = verbose,
                                trackB = trackB,
                                I_inv = I_inv,
                                Dy = Dy,
                                return_both_score_pvals = return_both_score_pvals,
                                cluster = cluster)
      
      if (is.null(test_result)) {
        if (return_nullB) {
          nullB_list[[test_ind]] <- NA
        }
      } else {
        
        score_test_hyperparams[test_ind, ] <- 
          c(test_result$u, test_result$rho, test_result$tau, test_result$inner_maxit,
            test_result$gap, test_result$convergence)
        
        if (return_nullB) {
          null_B <- test_result$null_B
          for (k in 1:p) {
            null_B[k, ] <- null_B[k, ] - constraint_fn(null_B[k, ])
          }
          nullB_list[[test_ind]] <- null_B
        }
        
        which_row <- which((as.numeric(coefficients$k) == as.numeric(test_kj$k[test_ind]))&
                             (as.numeric(coefficients$j) == as.numeric(test_kj$j[test_ind])))
        
        if (!return_both_score_pvals) {
          coefficients[which_row ,c("pval","score_stat")] <-
            c(test_result$pval,test_result$score_stat)
        } else {
          coefficients[which_row ,c("score_pval_full_info","score_stat_full_info",
                                    "score_pval_null_info","score_stat_null_info")] <-
            c(test_result$pval,test_result$score_stat,
              test_result$pval_null_info,test_result$score_stat_null_info)
        }
        
        if (use_both_cov) {
          
          #adjustment factor from guo GEE paper (https://doi.org/10.1002/sim.2161)
          alt_score_stat <- get_score_stat(Y = Y_test,
                                           X_cup = X_cup,
                                           X = X,
                                           B = test_result$null_B,
                                           k_constr = test_kj$k[test_ind],
                                           j_constr = test_kj$j[test_ind],
                                           constraint_grad_fn = constraint_grad_fn,
                                           indexes_to_remove = (j_ref - 1)*p + 1:p,
                                           j_ref = j_ref,
                                           J = J,
                                           n = n,
                                           p = p, 
                                           I_inv=I_inv,
                                           Dy = just_wald_things$Dy,
                                           cluster = cluster)
          
          
          which_row <- which((as.numeric(coefficients$k) == as.numeric(test_kj$k[test_ind]))&
                               (as.numeric(coefficients$j) == as.numeric(test_kj$j[test_ind])))
          coefficients[which_row, c("score_fullcov_p")] <- pchisq(alt_score_stat,1,
                                                                  lower.tail = FALSE)
        }
      
      }
    }
  }
  
  if (!is.null(colnames(X))) {
    if (length(unique(colnames(X))) == ncol(X)) {
      k_to_covariates <- data.frame(k = 1:p,
                                    covariate = colnames(X))
      
      coefficients$covariate <- do.call(c,
                                        lapply(  coefficients$k,
                                                 function(d) k_to_covariates$covariate[k_to_covariates$k ==d]))
    } else {
      coefficients$covariate <- NA
    }
  } else {
    coefficients$covariate <- NA
  }
  
  if (!is.null(colnames(Y))) {
    j_to_categories <- data.frame(j = 1:J,
                                  category = colnames(Y))
    coefficients$category <- do.call(c,
                                     lapply(coefficients$j,
                                            function(d) j_to_categories$category[
                                              j_to_categories$j ==d
                                            ]))
  } else {
    coefficients$category <- NA
  }
  
  if (!compute_cis) {
    coef_df <-
          cbind(data.frame(covariate = coefficients$covariate,
                           category = coefficients$category,
                           category_num = coefficients$j),
                coefficients[ , c("estimate","lower","upper")])
  } else {
    coef_df <-
      cbind(data.frame(covariate = coefficients$covariate,
                       category = coefficients$category,
                       category_num = coefficients$j),
            coefficients[ , c("estimate","se", "lower","upper")])
  }
  if (use_both_cov) {
    coef_df <- cbind(coef_df, coefficients[ , c("score_stat","pval","score_fullcov_p")])
  } else {
    if (return_both_score_pvals) {
      coef_df <- cbind(coef_df, coefficients[ , c("score_stat_full_info",
                                                  "score_pval_full_info",
                                                  "score_stat_null_info",
                                                  "score_pval_null_info")])
    } else {
      coef_df <- cbind(coef_df, coefficients[ , c("score_stat","pval")])
    }
  }
  if (return_wald_p) {
    coef_df$wald_p <- coefficients[, "wald_p"]
  }
  coefficients <- coef_df
  
  if (penalize) {
    if (!is.null(fitted_model)) {
      Y_augmented <- fitted_model$Y_augmented
    }
  } else {
    # set Y_augmented to NUll because without penalty there is no Y augmentation
    Y_augmented <- NULL
  }
  
  if (!is.null(fitted_model)) {
    if (penalize) {
      B <- fitted_model$B
    } else {
      B <- fitted_model
    }
  }
  
  if (penalize) {
    z_hat <- log(rowSums(Y_augmented)) - log(rowSums(exp(X %*% B)))
  } else {
    z_hat <- log(rowSums(Y)) - log(rowSums(exp(X %*% B)))
  }
  
  if (is.null(just_wald_things)) {
    I <- NULL
    Dy <- NULL
  } else {
    I <- just_wald_things$I
    Dy <- just_wald_things$Dy
  }
  
  if (is.null(colnames(B))) {
    colnames(B) <- colnames(Y)
  }
  
  if (is.null(rownames(B))) {
    rownames(B) <- colnames(X)
  }
  
  if (!is.null(zero_comparison_res)) {
    coefficients <- dplyr::full_join(coefficients, zero_comparison_res, 
                                     by = c("covariate", "category")) 
    coefficients$zero_comparison[is.na(coefficients$zero_comparison)] <- FALSE
    
    if (remove_zero_comparison_pvals == TRUE | is.numeric(remove_zero_comparison_pvals)) {
      pval_cols <- which(grepl("pval", names(coefficients)))
      for (col in pval_cols) {
        if (remove_zero_comparison_pvals == TRUE) {
          coefficients[coefficients$zero_comparison, col] <- NA
        } else {
          ind <- ifelse(is.na(coefficients[, col]), FALSE, 
                        coefficients[, col] <= remove_zero_comparison_pvals & coefficients$zero_comparison)
          coefficients[ind, col] <- NA
        }
      }
    }
  }
  
  results <- list("call" = call,
                  "coef" = coefficients,
                  "B" = B,
                  "penalized" = penalize,
                  "Y_augmented" = Y_augmented,
                  "z_hat" = z_hat,
                  "I" = I,
                  "Dy" = Dy,
                  "cluster" = cluster)
  
  if (refit & penalize) {
    results$estimation_converged <- converged_estimates
  }
  if (run_score_tests & return_nullB) {
    results$null_B <- nullB_list
  }
  if (run_score_tests) {
    results$score_test_hyperparams <- score_test_hyperparams
    if (sum(score_test_hyperparams$converged != "converged") > 0) {
      unconverged_test_kj <- test_kj[which(score_test_hyperparams$converged != "converged"), ]
      results$null_estimation_unconverged <- unconverged_test_kj
      warning("Optimization for estimation under the null for robust score tests failed to converge for some tests. See 'null_estimation_unconverged' within the returned emuFit object for which tests are affected by this.")
    }
  }
  
  return(structure(results, class = "emuFit"))
}


