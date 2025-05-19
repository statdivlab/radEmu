#' Fit radEmu model
#'
#' @param Y an n x J matrix or dataframe of nonnegative observations, or a phyloseq object containing an otu table and sample data.
#' @param X an n x p matrix or dataframe of covariates (optional)
#' @param formula a one-sided formula specifying the form of the mean model to be fit
#' @param data an n x p data frame containing variables given in \code{formula}
#' @param assay_name a string containing the desired assay name within a `TreeSummarizedExperiment` object.
#' This is only required if Y is a `TreeSummarizedExperiment` object, otherwise this argument does nothing
#' and can be ignored.
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
#' covariate k) of elements of B to construct hypothesis tests for. If you don't know
#' which indices k correspond to the covariate(s) that you would like to test, run the function
#' \code{radEmu::make_design_matrix()} in order to view the design matrix, and identify which
#' column of the design matrix corresponds to each covariate in your model. This argument is required when
#' running score tests.
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
#' @param match_row_names logical: Make sure rows on covariate data and response data correspond to 
#' the same sample by comparing row names and subsetting/reordering if necessary. 
#' @param constraint_fn function g defining a constraint on rows of B; g(B_k) = 0
#' for rows k = 1, ..., p of B. Default function is a smoothed median (minimizer of
#' pseudohuber loss). If a number is provided a single category constraint will be used
#' with the provided category as a reference category. This argument can either be a single constraint 
#' function to be used for all rows of B, or a list of length p of constraints to be used for each row of B.
#' @param constraint_grad_fn derivative of constraint_fn with respect to its
#' arguments (i.e., elements of a row of B). If \code{constraint_fn} is a list of constraint functions, then
#' this argument must also be a list. If \code{constraint_fn} is a single number, or a list that includes a 
#' single number, then the corresponding \code{constraint_grad_fn} can be set to \code{NULL}, and will be appropriately 
#' set within the function. 
#' @param constraint_param If pseudohuber centering is used (this is the default),
#' parameter controlling relative weighting of elements closer and further from center.
#' (Limit as \code{constraint_param} approaches infinity is the mean; as this parameter approaches zero,
#' the minimizer of the pseudo-Huber loss approaches the median.) If constraint function is not pseudohuber
#' centering (implemented in \code{radEmu::pseudohuber_median()}) then this argument will be ignored.
#' @param verbose provide updates as model is being fitted? Defaults to FALSE. If user sets verbose = TRUE,
#' then key messages about algorithm progress will be displayed. If user sets verbose = "development",
#' then key messages and technical messages about convergence will be displayed. Most users who want status
#' updates should set verbose = TRUE.
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
#' @param return_score_components logical: should components of score statistic be returned? Primarily used for debugging. Default is FALSE.
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
#' @param unobserved_taxon_error logical: should an error be thrown if Y includes taxa that have 0 counts for all samples? Default is TRUE.
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
#' @examples
#' # data frame example
#' data(wirbel_sample_small)
#' data(wirbel_otu_small)
#' emuRes <- emuFit(formula = ~ Group, data = wirbel_sample_small, Y = wirbel_otu_small,
#'                  test_kj = data.frame(k = 2, j = 1), tolerance = 0.01, 
#'                  constraint_tol = 0.01, B_null_tol = 0.01) 
#'  # here we set large tolerances for the example to run quickly, 
#'  # but we recommend smaller tolerances in practice
#'
#' # TreeSummarizedExperiment example (only run this if you have TreeSummarizedExperiment installed)
#' \dontrun{
#' library("TreeSummarizedExperiment")
#' example("TreeSummarizedExperiment")
#' assayNames(tse) <- "counts"
#' emuRes <- emuFit(Y = tse, formula = ~ condition, assay_name = "counts", 
#'                  test_kj = data.frame(k = 2, j = 1), tolerance = 0.01, constraint_tol = 0.01)
#'  # here we set large tolerances for the example to run quickly, 
#'  # but we recommend smaller tolerances in practice
#' }
#'
#' @export
#'
emuFit <- function(Y,
                   X = NULL,
                   formula = NULL,
                   data = NULL,
                   assay_name = NULL,
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
                   match_row_names = TRUE,
                   constraint_fn = pseudohuber_median,
                   constraint_grad_fn = dpseudohuber_median_dx,
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
                   return_score_components = FALSE,
                   return_both_score_pvals = FALSE,
                   remove_zero_comparison_pvals = 0.01,
                   unobserved_taxon_error = TRUE) {
  
  # Record call
  call <- match.call(expand.dots = FALSE)
  
  # run checks on arguments in function emuFit_check
  check_results <- emuFit_check(Y = Y,
                                X = X,
                                formula = formula,
                                data = data,
                                assay_name = assay_name,
                                cluster = cluster,
                                B_null_list = B_null_list,
                                test_kj = test_kj,
                                match_row_names = match_row_names,
                                verbose = verbose,
                                remove_zero_comparison_pvals = remove_zero_comparison_pvals,
                                unobserved_taxon_error = unobserved_taxon_error,
                                constraint_fn = constraint_fn,
                                constraint_grad_fn = constraint_grad_fn,
                                constraint_param = constraint_param,
                                run_score_tests = run_score_tests)
  
  Y <- check_results$Y
  X <- check_results$X
  cluster <- check_results$cluster
  B_null_list <- check_results$B_null_list
  test_kj <- check_results$test_kj
  constraint_fn <- check_results$constraint_fn
  constraint_grad_fn <- check_results$constraint_grad_fn
  constraint_param <- check_results$constraint_param
  
  n <- nrow(Y)
  J <- ncol(Y)
  p <- ncol(X)
  
  # check for zero-comparison parameters
  zero_comparison_res <- zero_comparison_check(X = X, Y = Y)
  
  X_cup <- X_cup_from_X(X,J)
  
  
  
  #choose ref taxon for fitting constrained models / performing wald and score tests
  j_ref <- get_j_ref(Y)
  
  ########################################
  ########## fit model if needed
  ########################################
  
  if (refit) {
    
    if (!is.null(fitted_model)) {
      B <- fitted_model$B
    }
    
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
                               verbose = (verbose == "development"),
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
                     constraint_fn = constraint_fn,
                     maxit = maxit,
                     max_stepsize = max_step,
                     tolerance = tolerance,
                     j_ref = j_ref,
                     verbose = (verbose == "development"))
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
    if (verbose %in% c(TRUE, "development")) {
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
                                    verbose = (verbose == "development"),
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
                                         converged = NA,
                                         niter = NA)
    
    if (return_score_components) {
      score_components <- vector(mode = "list", length = nrow(test_kj))
    }
    
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
    if (trackB) {
      trackB_list <- vector(mode = "list", length = nrow(test_kj))
    }
    for(test_ind in 1:nrow(test_kj)) {
      
      if (verbose %in% c(TRUE, "development")) {
        message(paste("Running score test ", test_ind, " of ", nrow(test_kj)," (row of B k = ", test_kj$k[test_ind], "; column of B j = ",
                    test_kj$j[test_ind],").",sep = ""))
        start <- proc.time()
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
                                verbose = (verbose == "development"),
                                trackB = trackB,
                                I_inv = I_inv,
                                Dy = Dy,
                                return_both_score_pvals = return_both_score_pvals,
                                cluster = cluster)
      
      if (return_score_components & !(is.null(test_result))) {
        score_components[[test_ind]] <- test_result$score_pieces
      }
      
      if (is.null(test_result)) {
        if (return_nullB) {
          nullB_list[[test_ind]] <- NA
        }
        if (trackB) {
          trackB_list[[test_ind]] <- NA
        }
      } else {
        
        score_test_hyperparams[test_ind, ] <- 
          c(test_result$u, test_result$rho, test_result$tau, test_result$inner_maxit,
            test_result$gap, test_result$convergence, test_result$niter)
        
        if (return_nullB) {
          null_B <- test_result$null_B
          for (k in 1:p) {
            null_B[k, ] <- null_B[k, ] - constraint_fn[[k]](null_B[k, ])
          }
          nullB_list[[test_ind]] <- null_B
        }
        
        if (trackB) {
          trackB_list[[test_ind]] <- test_result$Bs
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
                                           cluster = cluster)$score_stat
          
          
          which_row <- which((as.numeric(coefficients$k) == as.numeric(test_kj$k[test_ind]))&
                               (as.numeric(coefficients$j) == as.numeric(test_kj$j[test_ind])))
          coefficients[which_row, c("score_fullcov_p")] <- pchisq(alt_score_stat,1,
                                                                  lower.tail = FALSE)
        }
      
      }
      
      if (verbose %in% c(TRUE, "development")) {
        end <- proc.time() - start
        sec <- round(end[3]) 
        if (sec <= 300) {
          time <- paste0(sec, " seconds")
        } else if (sec <= 18000) {
          min <- round(sec / 60)
          time <- paste0(min, " minutes")
        } else {
          hour <- round(sec / (60^2))
          time <- paste0(hour, " hours")
        }
        message(paste("Score test ", test_ind, " of ", nrow(test_kj)," (row of B k = ", test_kj$k[test_ind], "; column of B j = ",
                    test_kj$j[test_ind],") has completed in approximately ", time, ".",sep = ""))
        start <- proc.time()
        
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
                  "cluster" = cluster,
                  "X" = X,
                  "constraint_fn" = constraint_fn,
                  "constraint_grad_fn" = constraint_grad_fn)
  
  if (is.null(Y_augmented)) {
    results$Y <- Y
  } 
  
  if (refit & penalize) {
    results$estimation_converged <- converged_estimates
  }
  if (run_score_tests) {
    if (return_nullB) {
      results$null_B <- nullB_list
    }
    if (trackB) {
      results$trackB_list <- trackB_list
    }
    results$score_test_hyperparams <- score_test_hyperparams
    if (sum(score_test_hyperparams$converged != "converged") > 0) {
      unconverged_test_kj <- test_kj[which(score_test_hyperparams$converged != "converged"), ]
      results$null_estimation_unconverged <- unconverged_test_kj
      warning("Optimization for estimation under the null for robust score tests failed to converge for some tests. See 'null_estimation_unconverged' within the returned emuFit object for which tests are affected by this.")
    }
  }
  if (run_score_tests & return_score_components) {
    results$score_components <- score_components
  }
  
  return(structure(results, class = "emuFit"))
}


