#' Generates the design matrix for a given \code{formula} and \code{data} object
#'
#' @param formula a one-sided formula specifying the form of the mean model to be fit
#' @param data an n x p data frame containing variables given in \code{formula} (required unless
#' \code{Y} is included as a \code{phyloseq} or \code{TreeSummarizedExperiment} object)
#' @param Y Optionally, a \code{phyloseq} or \code{TreeSummarizedExperiment} object 
#' containing an otu table and sample data.
#' @param assay_name Optionally, a string containing the desired assay name within a \code{TreeSummarizedExperiment} object.
#' This is only required if \code{data} is a \code{TreeSummarizedExperiment} object, otherwise this argument does nothing
#' and can be ignored.
#' 
#' @return returns design matrix \code{X}.
#' 
#' @importFrom stats model.matrix
#'
#' @export
#'
make_design_matrix <- function(Y = NULL,
                               formula,
                               data = NULL,
                               assay_name = NULL) {
  
  # check that Y or data is included
  if (is.null(Y) & is.null(data)) {
    stop("Data must be provided with the `data` or `Y` argument.")
  }
  
  # check if Y is a phyloseq object
  if ("phyloseq" %in% class(Y)) {
    if (requireNamespace("phyloseq", quietly = TRUE)) {
      data <- data.frame(phyloseq::sample_data(Y))
    } else {
      stop("You are trying to use a `phyloseq` data object or `phyloseq` helper function without having the `phyloseq` package installed. Please either install the package or use a standard data frame.")
    }
    
  # check if Y is a TreeSummarizedExperiment object
  } else if ("TreeSummarizedExperiment" %in% class(Y)) {
    if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      if (is.null(assay_name)) {
        stop("If Y is a `TreeSummarizedExperiment` object, make sure to include the assay_name and formula arguments.")
      }
      data <- as.data.frame(SummarizedExperiment::colData(Y))
    } else {
      stop("You are trying to use a `TreeSummarizedExperiment` data object or `TreeSummarizedExperiment` helper function without having the `SummarizedExperiment` package installed. Please either install the package or use a standard data frame.")
    }
  }
  
  X <- model.matrix(formula, data)
  return(X)
  
}