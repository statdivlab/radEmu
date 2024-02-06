#' Print function
#'
#' @param x Object of class \code{emuFit}
#' @param n The number of coefficient estimates to be printed (ordered by largest absolute value to smallest)
#' @param ... No optional arguments are accepted at this time.
#'
#' @return \code{NULL}. Displays printed model summary.
#'
#' @method print emuFit
#'
#' @export
print.emuFit <- function(x,
                         n = 20,
                         ...) {
  
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  
  cat("\nCoefficients estimates with the largest magnitudes:\n")
  n_or_nrow <- min(n, nrow(x$coef))
  sorted_ind <- order(abs(x$coef$estimate), decreasing = TRUE)[1:n_or_nrow]
  coefs_tab <- x$coef[sorted_ind,, drop = FALSE]
  cols_NA <- which(colMeans(is.na(x$coef)) == 1)
  if (length(cols_NA) > 0) {
    coefs_tab <- coefs_tab[, -(cols_NA)]
  }
  stats::printCoefmat(coefs_tab, na.print = "NA", ...)
  
  cat("\n")
  
  message("To obtain the entire coefficient table, use the command `emuFit_object$coef`.")
}
