#' Make lists of constraint functions and their gradients when using a reference taxon
#'
#' @param p The number of columns in the design matrix \code{X}. If you don't know the 
#' number of columns, you can find it with 
#' \code{ncol(radEmu::make_design_matrix(your_formula))}. \code{your_formula} should be 
#' the expression you give to \code{emuFit}'s \code{formula} argument.
#' @param j A single value or a vector of length \code{p - 1} where \code{p} is the number
#' of columns in the design matrix \code{X}. If a single value, \code{j} will be used as 
#' the reference category for all covariates. If a vector of values, \code{j[k]} will be
#' used as the reference category for the covariate in design matrix column \code{k + 1}.
#' 
#' @return A list with elements \code{constraints_list} and \code{constraints_grad_list}.
#' The \code{constraints_list} is a list of constraint functions for each column \code{p}
#' of the design matrix. By default, the constraint for the intercept is the pseudo Huber
#' median. The constraints for covariates are determined by reference categories given by
#' the argument \code{j}. The \code{constraints_grad_list} is a list of gradients of each
#' constraint function. 
#' 
#' @examples
#' # two columns in design matrix, reference taxon is taxon 5
#' list1 <- make_reference_constraints(p = 2, j = 5)
#' 
#' # four columns in design matrix, reference taxon for all covariates is taxon 2
#' list2 <- make_reference_constraints(p = 4, j = 2)
#' 
#' # four columns in design matrix, reference taxon for covariates 1 and 2 is taxon 3 and
#' # reference taxon for covariate 3 is taxon 4
#' list3 <- make_reference_constraints(p = 4, j = c(3, 3, 4))
#'
#' @export
#'
make_reference_constraints <- function(p, j) {
  if (length(j) > 1 & length(j) != p - 1) {
    stop("The argument `j` should be either a single value or a vector of length `p - 1`.")
  }
  constraint_list <- rep(list(NA), p)
  constraint_grad_list <- rep(list(NA), p)
  # set intercept constraint to be pseudo Huber median 
  constraint_list[[1]] <- (function(x) pseudohuber_center(x, d = 0.1))
  constraint_grad_list[[1]] <- (function(x) dpseudohuber_center_dx(x, d = 0.1))
  # set other constraints based on reference taxon 
  if (p > 2 & length(j) == 1) {
    j <- rep(j, p - 1)
  }
  for (k in 2:p) {
    constraint_cat <- j[k - 1]
    constraint_list[[k]] <- (function(constraint_cat) {
      force(constraint_cat)
      function(x) {x[constraint_cat]}})(constraint_cat)
    constraint_grad_list[[k]] <- (function(constraint_cat) {
      force(constraint_cat)
      function(x) {
        grad <- rep(0, length(x))
        grad[constraint_cat] <- 1
        return(grad)
      }})(constraint_cat)
  }
  # return results
  return(list(constraint_list = constraint_list, constraint_grad_list = constraint_grad_list))
}
