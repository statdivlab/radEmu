#' radEmu: Using Relative Abundance Data to Estimate of Multiplicative Differences in Mean Absolute Abundance 
#'
#' The radEmu package provides a method for estimating and testing multiplicative differences in mean absolute abundance, using relative abundance data.
#'
#' @docType package
#' @name radEmu
#' @keywords internal
"_PACKAGE"

## Silence R CMD check NOTE about undefined global variables
utils::globalVariables(c(
  "covariate", "estimate", "upper", "lower",
  "cat_small", "category", "constraint_diff",
  "it", "lik", "test_stat", "max_abs_B"
))