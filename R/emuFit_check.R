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
#' covariate k) of elements of B to construct hypothesis tests for. If \code{test_kj}
#' is not provided, all elements of B save the intercept row will be tested.
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
#' 
#' @return returns objects \code{Y}, \code{X}, \code{cluster}, and \code{B_null_list}, which may be modified by tests, and throw any useful
#' errors, warnings, or messages.
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
                         unobserved_taxon_error = TRUE) {
  
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
  
  # check for valid argument remove_zero_comparison_pvals 
  if (remove_zero_comparison_pvals != TRUE & remove_zero_comparison_pvals != FALSE) {
    if (!(is.numeric(remove_zero_comparison_pvals) & remove_zero_comparison_pvals <= 1 &
          remove_zero_comparison_pvals >= 0)) {
      stop("Please set `remove_zero_comparison_pvals` to either TRUE, FALSE, or a numeric value between 0 and 1.")
    }
  }
  
  return(list(Y = Y, X = X, cluster = cluster, B_null_list = B_null_list))
  
}