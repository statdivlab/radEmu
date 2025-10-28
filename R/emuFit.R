#' Fit radEmu model
#'
#' @param Y an n x J matrix or dataframe of nonnegative observations, or a phyloseq object containing an otu table and sample data.
#' @param X an n x p design matrix (either provide \code{X} or \code{data} and \code{formula})
#' @param formula a one-sided formula specifying the form of the mean model to be fit (used with \code{data})
#' @param data an n x p data frame containing variables given in \code{formula}
#' @param assay_name a string containing the desired assay name within a `TreeSummarizedExperiment` object.
#' This is only required if Y is a `TreeSummarizedExperiment` object, otherwise this argument can be ignored.
#' @param cluster a vector giving cluster membership for each row of Y to be used in computing 
#' GEE test statistics. Default is NULL, in which case rows of Y are treated as independent.
#' @param constraint_fn (Optional) User-provided constraint function, if default behavior of comparing log fold-difference
#' parameters to smoothed median over all categories is not desired. If a number is provided a single category constraint will be used
#' with the provided category as a reference category. This argument can either be a single constraint 
#' function to be used for all rows of B, or a list of length p of constraints to be used for each row of B.
#' @param constraint_grad_fn (Optional) User-provided derivative of constraint function, if default behavior of comparing log fold-difference
#' parameters to smoothed median over all categories is not desired. If \code{constraint_fn} is a list of constraint functions, then
#' this argument must also be a list. If \code{constraint_fn} is a single number, or a list that includes a 
#' single number, then the corresponding \code{constraint_grad_fn} can be set to \code{NULL}, and will be appropriately 
#' set within the function. 
#' @param constraint_param (Optional) If the smoothed median is used as a constraint (this is the default),
#' parameter controlling relative weighting of elements closer and further from center.
#' (Limit as \code{constraint_param} approaches infinity is the mean; as this parameter approaches zero,
#' the minimizer of the pseudo-Huber loss approaches the median.) If constraint function is not smoothed median
#'  (implemented in \code{radEmu::pseudohuber_median()}) then this argument will be ignored.
#' @param verbose provide updates as model is being fitted? Defaults to `FALSE`. If user sets `verbose = TRUE`,
#' then key messages about algorithm progress will be displayed. If user sets `verbose = "development"`,
#' then key messages and technical messages about convergence will be displayed. Most users who want status
#' updates should set `verbose = TRUE`.
#' @param match_row_names logical: If `TRUE`, make sure rows on covariate data and response data correspond to 
#' the same sample by comparing row names and subsetting/reordering if necessary. Default is `TRUE`.
#' @param unobserved_taxon_error logical: should an error be thrown if Y includes taxa that have 0 counts for all samples? Default is TRUE.
#' @param penalize logical: should Firth penalty be used? Default is `TRUE`. Used in estimation.
#' @param B starting value of coefficient matrix (p x J) for estimation. If not provided, B will be initiated as a zero matrix. 
#' Used in estimation.
#' @param fitted_model a fitted model produced by a separate call to emuFit; to be provided if score tests 
#' are to be run without refitting the full unrestricted model. Default is `NULL`.
#' @param refit logical: if `B` or `fitted_model` is provided, should full model be fit (`TRUE`) or
#' should fitting step be skipped (`FALSE`), e.g., if score tests are to be run on an already
#' fitted model. Default is `TRUE`.
#' @param tolerance tolerance for stopping criterion in estimation; once no element of B is updated by more than this value 
#' in a single step, we exit optimization. Defaults to `1e-3`. Used in estimation. 
#' @param maxit maximum number of outer iterations to perform before exiting optimization. Default is `1000`.
#' Used in estimation.
#' @param alpha nominal type 1 error level to be used to construct confidence intervals. Default is `0.05`
#' (corresponding to 95% confidence intervals)
#' @param return_wald_p logical: return p-values from Wald tests? Default is `FALSE`.
#' @param compute_cis logical: compute and return Wald CIs? Default is `TRUE`.
#' @param run_score_tests logical: perform robust score testing? Default is TRUE.
#' @param test_kj a data frame whose rows give coordinates (in category j and
#' covariate k) of elements of B to construct hypothesis tests for. If you don't know
#' which indices k correspond to the covariate(s) that you would like to test, run the function
#' \code{radEmu::make_design_matrix()} in order to view the design matrix, and identify which
#' column of the design matrix corresponds to each covariate in your model. This argument is required when
#' running score tests.
#' @param null_fit_alg Which null fitting algorithm to use for score tests: \code{"constraint_sandwich"} or 
#' \code{"augmented_lagrangian"}. Default and recommended approach is \code{"constraint_sandwich"}, unless \code{J < 20}.
#' @param B_null_list list of starting values of coefficient matrix (p x J) for null estimation for score testing. This should either 
#' be a list with the same length as \code{test_kj}. If you only want to provide starting values for some tests,
#' include the other elements of the list as \code{NULL}.
#' @param maxit_null maximum number of outer iterations to perform before exiting optimization. Default is `1000`.
#' Used in estimation under null hypothesis for score tests.
#' @param tol_lik tolerance for relative changes in likelihood for stopping criteria. Default is `1e-5`. Used in 
#' estimation under null hypothesis for score tests with "constraint_sandwich" algorithm. 
#' @param tol_test_stat tolerance for relative changes in test statistic for stopping criteria. Default is `0.01`. Used in
#' estimation under null hypothesis for score tests with "constraint_sandwich" algorithm. 
#' @param null_window window to use for stopping criteria (this many iterations where stopping criteria is met). Default is `5`.
#' Used in estimation under null hypothesis for score tests with "constraint_sandwich" algorithm.
#' @param null_diagnostic_plots logical: should diagnostic plots be made for estimation under the null hypothesis? Default is \code{FALSE}.
#' @param remove_zero_comparison_pvals Should score p-values be replaced with NA for zero-comparison parameters? These parameters occur 
#' for categorical covariates with three or more levels, and represent parameters that compare a covariate level to the reference level for
#' a category in which the comparison level and reference level both have 0 counts in all samples. These parameters can have misleadingly 
#' small p-values and are not thought to have scientifically interesting signals. We recommend removing them before analyzing data further. 
#' If TRUE, all zero-comparison parameter p-values will be set to NA. If FALSE no zero-comparison parameter p-values will be set to NA.
#' If a value between 0 and 1, all zero-comparison p-values below the value will be set to NA. 
#' Default is \code{0.01}. 
#' @param control A list of control parameters, to have more control over estimation and hypothesis testing. See \code{control_fn} for details.
#' @param ... Additional arguments. Arguments matching the names of \code{control_fn()} options are forwarded to that function and override
#' defaults. Unknown arguments are ignored with a warning.
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
#'                  test_kj = data.frame(k = 2, j = 1), tolerance = 0.01) 
#'  # here we set large tolerances for the example to run quickly, 
#'  # but we recommend smaller tolerances in practice
#'
#' # TreeSummarizedExperiment example (only run this if you have TreeSummarizedExperiment installed)
#' \dontrun{
#' library("TreeSummarizedExperiment")
#' example("TreeSummarizedExperiment")
#' assayNames(tse) <- "counts"
#' emuRes <- emuFit(Y = tse, formula = ~ condition, assay_name = "counts", 
#'                  test_kj = data.frame(k = 2, j = 1), tolerance = 0.01)
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
                   constraint_fn = pseudohuber_median,
                   constraint_grad_fn = dpseudohuber_median_dx,
                   constraint_param = 0.1,
                   verbose = FALSE,
                   match_row_names = TRUE,
                   unobserved_taxon_error = TRUE,
                   penalize = TRUE,
                   B = NULL,
                   fitted_model = NULL,
                   refit = TRUE,
                   tolerance = 1e-4,
                   maxit = 1000,
                   alpha = 0.05,
                   return_wald_p = FALSE,
                   compute_cis = TRUE,
                   run_score_tests = TRUE,
                   test_kj = NULL,
                   null_fit_alg = "constraint_sandwich",
                   B_null_list = NULL,
                   maxit_null = 1000,
                   tol_lik = 1e-5,
                   tol_test_stat = 0.01,
                   null_window = 5,
                   null_diagnostic_plots = FALSE, 
                   remove_zero_comparison_pvals = 0.01,
                   control = NULL,
                   ...) {
  
  # Record call
  call <- match.call(expand.dots = FALSE)
  
  # capture ...
  dots <- list(...)
  
  # get possible control args from formals(control_fn)
  control_args <- names(formals(control_fn))
  control_args <- setdiff(control_args, "control") 
  
  # split into control_dots vs unknown
  control_dots <- dots[names(dots) %in% control_args]
  unknown_dots <- dots[!names(dots) %in% control_args]
  
  if (length(unknown_dots) > 0) {
    warning("Unknown arguments in emuFit(...): ", paste(names(unknown_dots), collapse = ", "))
  }
  
  # merge with user-supplied control list
  if (is.null(control)) {
    control <- control_fn(control_dots)
  } else {
    control <- control_fn(utils::modifyList(control, control_dots))
  }
  
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
                                run_score_tests = run_score_tests,
                                null_fit_alg = null_fit_alg)
  
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
  if (J < 20) {
    null_fit_alg <- "augmented_lagrangian"
  }
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
                               max_step = control$max_step,
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
                     max_stepsize = control$max_step,
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
    if (control$use_fullmodel_cov) {
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
  if (control$use_both_cov) {
    coefficients$score_fullcov_p <- NA
  }
  
  if (control$use_fullmodel_info) {
    indexes_to_remove = (j_ref - 1)*p + 1:p
    I_inv <- Matrix::solve(just_wald_things$I[-indexes_to_remove,-indexes_to_remove],  method = "cholmod_solve")
  } else {
    if(control$use_both_cov){
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
    
    if (control$return_score_components) {
      score_components <- vector(mode = "list", length = nrow(test_kj))
    }
    
    if (control$return_both_score_pvals) {
      colnames(coefficients)[colnames(coefficients) == "pval"] <-
        "score_pval_full_info"
      colnames(coefficients)[colnames(coefficients) == "score_stat"] <-
        "score_stat_full_info"
      coefficients$score_pval_null_info <- NA
      coefficients$score_stat_null_info <- NA
      
      if (!control$use_fullmodel_info) {
        stop("If control$return_both_score_pvals = TRUE, control$use_fullmodel_info must be TRUE as well.")
      }
    }
    
    if (control$return_nullB) {
      nullB_list <- vector(mode = "list", length = nrow(test_kj))
    }
    if (control$trackB) {
      trackB_list <- vector(mode = "list", length = nrow(test_kj))
    }
    if (null_diagnostic_plots) {
      null_plots <- vector(mode = "list", length = nrow(test_kj))
    }
    
    if (null_fit_alg == "constraint_sandwich") {
      k_list <- list()
      for (k in unique(test_kj$k)) {
        # set as other by default
        constraint_type <- "other"
        
        # check if it is a single category constraint
        v1 <- 1:J
        v2 <- c(J, 1:(J - 1))
        r <- constraint_fn[[k]](v1)
        s <- constraint_fn[[k]](v2)
        r_grad <- constraint_grad_fn[[k]](v1)
        s_grad <- constraint_grad_fn[[k]](v2)
        if ((r == 1 && s == J) | r > 1 && s == (r - 1)) {
          expected_r_grad <- rep(0, J)
          expected_r_grad[v1 == r] <- 1
          if (isTRUE(all.equal(r_grad, expected_r_grad)) &&
              isTRUE(all.equal(r_grad, s_grad))) {
            constraint_type <- "scc"
          }
        } 
        
        # check if it is symmetric constraint 
        # first check mean 
        v3 <- rnorm(J)
        if (constraint_fn[[k]](v3) == mean(v3)) {
          constraint_type <- "symmetric"
        }
        # then check pseudo-Huber median 
        if (any(grepl("pseudohuber_median", deparse(body(constraint_fn[[k]]))))) {
          constraint_type <- "symmetric"
        }
        
        # check if it is symmetric constraint over a subset 
        ref_set <- try(get("reference_set", envir = environment(constraint_fn[[k]]), inherits = TRUE))
        uses_ref <- any(grepl("reference_set", deparse(body(constraint_fn[[k]]))))
        if (!inherits(ref_set, "try-error") & uses_ref) {
          if (constraint_fn[[k]](v3) == mean(v3[ref_set[[k]]])) {
            constraint_type <- "symmetric_subset"
          }
          if (any(grepl("pseudohuber_median", deparse(body(constraint_fn[[k]]))))) {
            constraint_type <- "symmetric_subset"
          }
        }

        k_list[[k]] <- constraint_type
      }
    } else {
      k_list <- as.list(rep("other", max(test_kj$k)))
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
                                rho_init = control$rho_init,
                                tau = control$tau,
                                kappa = control$kappa,
                                B_tol = control$B_null_tol,
                                inner_tol = control$inner_tol,
                                constraint_tol = control$constraint_tol,
                                j_ref = j_ref,
                                c1 = control$c1,
                                maxit = maxit_null,
                                inner_maxit = control$inner_maxit,
                                ntries = control$ntries,
                                verbose = (verbose == "development"),
                                trackB = control$trackB,
                                I_inv = I_inv,
                                Dy = Dy,
                                return_both_score_pvals = control$return_both_score_pvals,
                                cluster = cluster,
                                null_fit_constraint = k_list[[test_kj$k[test_ind]]],
                                null_diagnostic_plots = null_diagnostic_plots,
                                ignore_stop = control$ignore_stop,
                                tol_lik = tol_lik,
                                tol_test_stat = tol_test_stat,
                                null_window = null_window)
      
      if (is.null(test_result)) {
        if (control$return_nullB) {
          nullB_list[[test_ind]] <- NA
        }
        if (control$trackB) {
          trackB_list[[test_ind]] <- NA
        }
      } else {
        
        if (control$return_score_components) {
          score_components[[test_ind]] <- test_result$score_pieces
        }
        
        if (null_diagnostic_plots) {
          null_plots[[test_ind]] <- test_result$diagnostics
        }
        
        cols_to_update <- intersect(names(score_test_hyperparams), names(test_result))
        score_test_hyperparams[test_ind, cols_to_update] <- test_result[cols_to_update]
        
        # score_test_hyperparams[test_ind, ] <- 
        #   c(test_result$u, test_result$rho, test_result$tau, test_result$inner_maxit,
        #     test_result$gap, test_result$convergence, test_result$niter)
        
        if (control$return_nullB) {
          null_B <- test_result$null_B
          for (k in 1:p) {
            null_B[k, ] <- null_B[k, ] - constraint_fn[[k]](null_B[k, ])
          }
          nullB_list[[test_ind]] <- null_B
        }
        
        if (control$trackB) {
          trackB_list[[test_ind]] <- test_result$Bs
        }
        
        which_row <- which((as.numeric(coefficients$k) == as.numeric(test_kj$k[test_ind]))&
                             (as.numeric(coefficients$j) == as.numeric(test_kj$j[test_ind])))
        
        if (!control$return_both_score_pvals) {
          coefficients[which_row ,c("pval","score_stat")] <-
            c(test_result$pval,test_result$score_stat)
        } else {
          coefficients[which_row ,c("score_pval_full_info","score_stat_full_info",
                                    "score_pval_null_info","score_stat_null_info")] <-
            c(test_result$pval,test_result$score_stat,
              test_result$pval_null_info,test_result$score_stat_null_info)
        }
        
        if (control$use_both_cov) {
          
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
  if (control$use_both_cov) {
    coef_df <- cbind(coef_df, coefficients[ , c("score_stat","pval","score_fullcov_p")])
  } else {
    if (control$return_both_score_pvals) {
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
  if (run_score_tests) {
    if (control$return_nullB) {
      results$null_B <- nullB_list
    }
    if (control$trackB) {
      results$trackB_list <- trackB_list
    }
    results$score_test_hyperparams <- score_test_hyperparams
    if (sum(score_test_hyperparams$converged != "converged") > 0) {
      unconverged_test_kj <- test_kj[which(score_test_hyperparams$converged != "converged"), ]
      results$null_estimation_unconverged <- unconverged_test_kj
      warning("Optimization for estimation under the null for robust score tests failed to converge for some tests. See 'null_estimation_unconverged' within the returned emuFit object for which tests are affected by this.")
    }
    
    if (control$return_score_components) {
      results$score_components <- score_components
    }
    
    if (null_diagnostic_plots) {
      results$null_diagnostic_plots <- null_plots
    }
  }
  
  return(structure(results, class = "emuFit"))
}


