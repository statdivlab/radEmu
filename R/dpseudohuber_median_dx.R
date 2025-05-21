#' Calculate the derivative to the pseudo-Huber smoothed median
#'
#' @param x A vector to calculate the derivative of the pseudo-Huber smoothed median for. 
#' @param d Smoothing parameter, by default set to \code{0.1}. As \code{d} approaches 
#' \code{0} the pseudo-Huber median function approaches the median and as \code{d} approaches 
#' infinity this function approaches the mean. 
#' @param na.rm Passed to \code{pseudohuber_median}, default is FALSE, if FALSE then when 
#' \code{x} includes at least one NA value then NA is returned, if TRUE then when \code{x} 
#' includes at least one NA value then that value is removed and the pseudo-Huber median 
#' is computed without it.
#'
#' @return The derivative of the calculated pseudo-Huber smoothed median over \code{x} with smoothing 
#' parameter \code{d}. 
#' 
#' @examples
#' dpseudohuber_median_dx(x = rnorm(10), d = 0.1)
#'
#' @aliases dpsuedohuber_median_dx
#'
#' @export
#'
dpseudohuber_median_dx <- function(x,
                                   d = 0.1,
                                   na.rm = FALSE){
  ps_center <- pseudohuber_median(x,d, na.rm = na.rm)
  scaled_sq <- ((x - ps_center)/d)^2
  #derivation of why the following is the derivative of 
  #the pseudohuber centering is given in supplement of 
  #Willis & Clausen (2024)
  w <- sqrt(1/(1 + scaled_sq))
  w3 <- w^3
  return((w3/sum(w3, na.rm = na.rm)))
}

#' @rdname dpseudohuber_median_dx
#' @export
dpsuedohuber_median_dx <- dpseudohuber_median_dx 
