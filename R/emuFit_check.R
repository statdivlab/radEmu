#' Runs checks for appropriate arguments before running \code{emuFit()}
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
#' @param B_null_list list of starting values of coefficient matrix (p x J) for null estimation. This should either 
#' be a list with the same length as \code{test_kj}. If you only want to provide starting values for some tests,
#' include the other elements of the list as \code{NULL}.
#' @param test_kj a data frame whose rows give coordinates (in category j and
#' covariate k) of elements of B to construct hypothesis tests for. `k` could also be the name of a covariate 
#' included in `X` or `data`. If you don't know which coordinates k correspond to the covariate(s) that you 
#' would like to test, run the function \code{radEmu::make_design_matrix()} in order to view the design matrix,
#' and identify which column of the design matrix corresponds to each covariate in your model. This argument 
#' is required when running score tests.
#' @param match_row_names logical: Make sure rows on covariate data and response data correspond to 
#' the same sample by comparing row names and subsetting/reordering if necessary. 
#' @param verbose provide updates as model is being fitted? Defaults to FALSE. If user sets verbose = TRUE,
#' then key messages about algorithm progress will be displayed. If user sets verbose = "development",
#' then key messages and technical messages about convergence will be displayed. Most users who want status
#' updates should set verbose = TRUE.
#' @param remove_zero_comparison_pvals Should score p-values be replaced with NA for zero-comparison parameters? These parameters occur 
#' for categorical covariates with three or more levels, and represent parameters that compare a covariate level to the reference level for
#' a category in which the comparison level and reference level both have 0 counts in all samples. These parameters can have misleadingly 
#' small p-values and are not thought to have scientifically interesting signals. We recommend removing them before analyzing data further. 
#' If TRUE, all zero-comparison parameter p-values will be set to NA. If FALSE no zero-comparison parameter p-values will be set to NA.
#' If a value between 0 and 1, all zero-comparison p-values below the value will be set to NA. 
#' Default is \code{0.01}. 
#' @param unobserved_taxon_error logical: should an error be thrown if Y includes taxa that have 0 counts for all samples? Default is TRUE.
#' @param constraint_fn function g defining a constraint on rows of B; g(B_k) = 0
#' for rows k = 1, ..., p of B. Default function is a smoothed median (minimizer of
#' pseudohuber loss). If a number is provided a single category constraint will be used
#' with the provided category as a reference category. This argument can either be a single constraint 
#' function to be used for all rows of B, or a list of length p of constraints to be used for each row of B.
#' @param constraint_grad_fn derivative of constraint_fn with respect to its
#' arguments (i.e., elements of a row of B). If \code{constraint_fn} is a list of constraint functions, then
#' this argument must also be a list.
#' @param constraint_param If pseudohuber centering is used (this is the default),
#' parameter controlling relative weighting of elements closer and further from center.
#' (Limit as \code{constraint_param} approaches infinity is the mean; as this parameter approaches zero,
#' the minimizer of the pseudo-Huber loss approaches the median.)
#' @param run_score_tests logical: perform robust score testing? 
#' @param null_fit_alg Which null fitting algorithm to use \code{"fisher_scoring"} or \code{"augmented_lagrangian"}.
#' Default and recommended approach is \code{"fisher_scoring"}.
#' 
#' @return returns objects \code{Y}, \code{X}, \code{cluster}, and \code{B_null_list}, which may be modified by tests, and throw any useful
#' errors, warnings, or messages.
#' 
#' @importFrom stats model.matrix
#'
emuFit_check <- function(Y,
                         X = NULL,
                         formula = NULL,
                         data = NULL,
                         assay_name = NULL,
                         cluster = NULL,
                         B_null_list = NULL,
                         test_kj = NULL,
                         match_row_names = TRUE,
                         verbose = FALSE,
                         remove_zero_comparison_pvals = 0.01,
                         unobserved_taxon_error = TRUE,
                         constraint_fn,
                         constraint_grad_fn,
                         constraint_param,
                         run_score_tests = TRUE,
                         null_fit_alg = "constraint_sandwich") {
  
  # confirm that input to verbose is valid
  if (!(verbose %in% c(FALSE, TRUE, "development"))) {
    stop('The argument "verbose" must be set to one of TRUE, FALSE, or "development".')
  }
  
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
    
    # check if Y is a TreeSummarizedExperiment object
  } else if ("TreeSummarizedExperiment" %in% class(Y)) {
    if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      if (is.null(assay_name) | is.null(formula)) {
        stop("If Y is a `TreeSummarizedExperiment` object, make sure to include the assay_name and formula arguments.")
      }
      data <- as.data.frame(SummarizedExperiment::colData(Y))
      X <- model.matrix(formula, data)
      Y <- as.data.frame(t(SummarizedExperiment::assay(Y, assay_name)))
    } else {
      stop("You are trying to use a `TreeSummarizedExperiment` data object or `TreeSummarizedExperiment` helper function without having the `SummarizedExperiment` package installed. Please either install the package or use a standard data frame.")
    }
  }
  
  # convert Y from a data.frame object to a matrix, even if it was extracted directly from `phyloseq` or `TreeSummarizedExperiment`
  if ("data.frame" %in% class(Y)) {
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
  
  # check for NA's in X or that have been removed going from data to X
  if (is.null(data)) {
    if (sum(is.na(X)) > 0) {
      stop("At least one value in your design matrix `X` is NA. Please remove any observations with missing values (NAs) or impute these missing values according to your analysis plan.")
    }
  } else {
    if (nrow(X) < nrow(data)) {
      warning("Your data includes at least one NA value, and the `model.matrix()` function is automatically dropping all observations with any missing values. If you want to avoid this default behavior, please manually remove observations from your data with missing values or impute these missing values according to your analysis plan.")
    }
  }
  
  # check that if X and Y match in the row names
  if (is.null(rownames(X)) || is.null(rownames(Y))){
    if (nrow(X) == nrow(Y)){
      if(match_row_names){
        if(is.null(rownames(X))){
          message("Row names are missing from the covariate matrix X. We will assume the rows are in the same order as in the response matrix Y. You are responsible for ensuring the order of your observations is the same in both matrices.")
        } else {
          message("Row names are missing from the response matrix Y. We will assume the rows are in the same order as in the covariate matrix X. You are responsible for ensuring the order of your observations is the same in both matrices.")
        }
      }
    } else {
      if(is.null(rownames(X))){
        stop("Row names are missing from the covariate matrix X, and the number of rows does not match the number of rows in the response matrix Y. Please resolve this issue before refitting the model.")
      } else {
        stop("Row names are missing from the response matrix Y, and the number of rows does not match the number of rows in the covariate matrix X. Please resolve this issue before refitting the model.")
      }
    }
  } else{
    if(match_row_names){
      names_X <- rownames(X)
      names_Y <- rownames(Y)
      
      #Checking if any row names are duplicated
      if (any(duplicated(names_X))) stop("Covariate matrix X has duplicated row names. Please ensure all row names are unique before refitting the model.")
      if (any(duplicated(names_Y))) stop("Response matrix Y has duplicated row names. Please ensure all row names are unique before refitting the model.")
      
      # Find common row names
      common_names <- intersect(names_X, names_Y)
      
      if (length(common_names) < length(names_X) || length(common_names) < length(names_Y)) {
        warning(sprintf("According to the rownames, there are observations that are missing either in the covariate matrix (X) and/or the response matrix (Y). We will subset to common rows only, resulting in %d samples.", length(common_names))) 
        
        X <- X[common_names, , drop = FALSE]
        Y <- Y[common_names, , drop = FALSE]
        
      } else if(all.equal(rownames(Y), rownames(X)) != TRUE){
        message("There is a different row ordering between the covariate matrix (X) and the response matrix (Y). Covariate data will be reordered to match response data.")
        X <- X[rownames(Y), , drop = FALSE]
      }
    } else {
      if(nrow(X) != nrow(Y)){
        stop("The number of rows does not match between the covariate matrix (X) and the response matrix (Y), and subsetting/matching by row name has been disabled. Please resolve this issue before refitting the model.")
      }
    }
  }
  
  p <- ncol(X)
  J <- ncol(Y)
  if (is.null(colnames(X))) {
    if (p > 1) {
      colnames(X) <- c("Intercept", paste0("covariate_", 1:(p - 1)))
    } else {
      colnames(X) <- "Intercept"
    }
  }
  if (is.null(colnames(Y))) {
    colnames(Y) <- paste0("category_", 1:J)
  }
  
  if (min(rowSums(Y))==0) {
    stop("Some rows of Y consist entirely of zeroes, meaning that some samples
have no observations. These samples must be excluded before fitting model.")
  }
  
  if (unobserved_taxon_error) {
    if (min(colSums(Y)) == 0) {
      stop("Some columns of Y consist entirely of zeroes, meaning that some categories have zero counts for all samples. These
         categories must be excluded before fitting the model.")
    }
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
  
  # check that test kj values are numbers, if they are strings that convert if possible to numbers
  if (!is.null(test_kj)) {
    
    if (!is.numeric(test_kj$k)) {
      cn <- colnames(X)
      unique_k <- unique(test_kj$k)
      
      k_map <- lapply(unique_k, function(k_j) {
        
        num <- suppressWarnings(as.numeric(k_j))
        if (!is.na(num)) return(num)
        
        exact_idx <- which(cn == k_j)
        if (length(exact_idx) == 1) return(exact_idx)
        
        partial_idx <- which(startsWith(cn, k_j))
        if (length(partial_idx) > 0) return(partial_idx)
        
        stop("Make sure that the values of `k` in `test_kj` are numeric or correspond to column names of the `X` matrix.")
        
      })
      
      names(k_map) <- unique_k
      expanded_list <- vector("list", length = nrow(test_kj) * ncol(X))
      idx <- 1
      
      for (row in seq_len(nrow(test_kj))) {
        k_j <- test_kj$k[row]
        matches <- k_map[[as.character(k_j)]]
        test_kj$k[row] <- matches[1]
        if (length(matches) > 1) {
          for (m in 2:length(matches)) {
            tmp <- test_kj[row, ]
            tmp$k <- matches[m]
            expanded_list[[idx]] <- tmp
            idx <- idx + 1
          }
        }
      }
      
      test_kj_expanded <- do.call(rbind, expanded_list[seq_len(idx - 1)])
      test_kj <- rbind(test_kj, test_kj_expanded)
      test_kj$k <- as.numeric(test_kj$k)
    }
    # if (!(is.numeric(test_kj$k))) {
    #   if (sum(!(test_kj$k %in% colnames(X))) == 0) {
    #     test_kj$k <- as.vector(sapply(test_kj$k, function(x) {which(colnames(X) == x)}))
    #   } else {
    #     stop("Make sure that the values of `k` in `test_kj` are numeric or correspond to column names of the `X` matrix.")
    #   }
    # }
    if (!(is.numeric(test_kj$j))) {
      if (sum(!(test_kj$j %in% colnames(Y))) == 0) {
        test_kj$j <- as.vector(sapply(test_kj$j, function(x) {which(colnames(Y) == x)}))
      } else {
        stop("Make sure that the values of `j` in `test_kj` are numeric or correspond to column names of the `Y` matrix.")
      }
    }
  }
  
  # check that test_kj is not null if running score tests
  if (run_score_tests) {
    if (is.null(test_kj)) {
      stop("When `run_score_tests = TRUE`, you must provide a matrix `test_kj` to determine which parameters you want to test. If you don't know which indices k correspond to the covariate(s) that you would like to test, run the function `radEmu::make_design_matrix()` in order to view the design matrix, and identify which column of the design matrix corresponds to each covariate in your model. If you don't know which indices j correspond to categories (taxa) that you want to test, you can look at the columns and column names of your `Y` matrix.")
    }
  }
  
  # check for valid argument remove_zero_comparison_pvals 
  if (remove_zero_comparison_pvals != TRUE & remove_zero_comparison_pvals != FALSE) {
    if (!(is.numeric(remove_zero_comparison_pvals) & remove_zero_comparison_pvals <= 1 &
          remove_zero_comparison_pvals >= 0)) {
      stop("Please set `remove_zero_comparison_pvals` to either TRUE, FALSE, or a numeric value between 0 and 1.")
    }
  }
  
  # check for valid constraints
  if (length(constraint_fn) != length(constraint_grad_fn)) {
    stop("The arguments `constraint_fn` and `constraint_grad_fn` have different lengths. These arguments should either both be single functions, or both be lists of p functions, for p rows of the B matrix.")
  }
  if (length(constraint_fn) == 1) {
    if (inherits(constraint_fn, "list")) {
      constraint_fn <- constraint_fn[[1]]
    }
    if (inherits(constraint_grad_fn, "list")) {
      constraint_grad_fn <- constraint_grad_fn[[1]]
    }
    constraint_fn <- rep(list(constraint_fn), p)
    constraint_grad_fn <- rep(list(constraint_grad_fn), p)
  }

  if (length(constraint_fn) != p) {
    stop("The arguments `constraint_fn` and `constraint_grad_fn` have the wrong lengths. These arguments should either both be single functions, or both be lists of p functions, for p rows of the B matrix.")
  }
  
  for (k in 1:p) {
    if (length(constraint_fn[[k]]) == 1 & is.numeric(constraint_fn[[k]])) {
      constraint_cat <- constraint_fn[[k]]
      constraint_fn[[k]] <- (function(constraint_cat) {
        force(constraint_cat)
        function(x) {x[constraint_cat]}})(constraint_cat)
      constraint_grad_fn[[k]] <- (function(constraint_cat) {
        force(constraint_cat)
        function(x) {
        grad <- rep(0, length(x))
        grad[constraint_cat] <- 1
        return(grad)
      }})(constraint_cat)
    }
    
    if (is.logical(all.equal(constraint_fn[[k]], pseudohuber_median))) {
      if (all.equal(constraint_fn[[k]], pseudohuber_median)) {
        if (verbose %in% c(TRUE, "development")) message("Centering row ", k, " of B with pseudo-Huber smoothed median with smoothing parameter ", constraint_param, ".")
        
        stopifnot(!is.na(constraint_param))
        
        constraint_fn[[k]] <- (function(x) pseudohuber_median(x, d = constraint_param))
        constraint_grad_fn[[k]] <- (function(x) dpseudohuber_median_dx(x, d = constraint_param))
        
      } 
    } else {
      
      # commenting out this warning. Adding to the documentation to reflect that this parameter
      # will only be used when centering with pseudohuber center function
      
#       if (!is.na(constraint_param)) {
#         
#         warning("Argument constraint_param is currently only supported for centering with
# pseudohuber_center() function; constraint_param input is otherwise ignored. Please directly
# feed your choice of constraint function, including any necessary parameters, to constraint_fn argument
# and the corresponding gradient function to constraint_grad_fn.")
#       }
    } 
  }
 
  # check that an appropriate null_fit_alg is given
  if (!(null_fit_alg %in% c("constraint_sandwich", "augmented_lagrangian"))) {
    stop("The two options for `null_fit_alg` are 'constraint_sandwich' or 'augmented_lagrangian'.")
  }
  
  return(list(Y = Y, X = X, cluster = cluster, B_null_list = B_null_list, test_kj = test_kj,
              constraint_fn = constraint_fn, constraint_grad_fn = constraint_grad_fn, 
              constraint_param = constraint_param))
  
}
