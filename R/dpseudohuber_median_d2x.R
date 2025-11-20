#' Calculate the second derivative of the pseudo-Huber smoothed median
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
#' @return The second derivative of the calculated pseudo-Huber smoothed median over \code{x} with smoothing 
#' parameter \code{d}. 
#' 
#' @note algebra derived by ChatGPT based on our supp
#' 
#' @examples
#' dpseudohuber_median_dx(x = rnorm(10), d = 0.1)
#'
#' @aliases dpsuedohuber_median_dx
#'
#' @export
#'
dpseudohuber_median_d2x <- function(x,
                                    d = 0.1,
                                    na.rm = FALSE) {
  # pseudo-Huber center
  ps_center <- pseudohuber_median(x, d, na.rm = na.rm)
  
  # basic quantities
  y  <- (x - ps_center) / d
  w  <- 1 / sqrt(1 + y^2)          # w_j
  w3 <- w^3                        # a_j = w_j^3
  S  <- sum(w3, na.rm = na.rm)     # S = sum_j w_j^3
  
  # t_j = (x_j - c*)/d^2 * w_j^5
  v  <- (x - ps_center) / d^2
  w5 <- w^5
  t  <- v * w5                     # t_j
  T  <- sum(t, na.rm = na.rm)      # T = sum_j t_j
  
  n <- length(x)
  H <- matrix(NA_real_, n, n)
  
  # Hessian entries:
  # H_{pq} = -3/S * [ t_p (delta_{pq} - w3_q/S) - (t_q - w3_q*T/S) * w3_p/S ]
  for (p in seq_len(n)) {
    a_p <- w3[p]
    t_p <- t[p]
    for (q in seq_len(n)) {
      a_q <- w3[q]
      t_q <- t[q]
      delta_pq <- as.numeric(p == q)
      
      H[p, q] <- -(3 / S) * (
        t_p * (delta_pq - a_q / S) -
          (t_q - a_q * T / S) * a_p / S
      )
    }
  }
  
  H
}