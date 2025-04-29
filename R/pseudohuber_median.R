#' Calculate the pseudo-Huber smoothed median
#'
#' Calculate the pseudo-Huber smoothed median using a quadratic approximation to the 
#' pseudo-Huber criterion detailed in supplement to Clausen et al. (2024).
#'
#' @param x A vector to calculate the pseudo-Huber smoothed median for. 
#' @param d Smoothing parameter, by default set to \code{0.1}. As \code{d} approaches 
#' \code{0} this function approaches the median and as \code{d} approaches infinity
#' this function approaches the mean. 
#' @param tolerance Tolerance used to determine convergence in the algorithm used 
#' to calculate this function value. 
#'
#' @return The calculated pseudo-Huber smoothed median over \code{x} with smoothing 
#' parameter \code{d}. 
#' 
#' @examples
#' pseudohuber_median(x = rnorm(10), d = 0.1)
#' 
#' @aliases psuedohuber_median
#' 
#' @export
#'
pseudohuber_median <- function(x,
                               d = 0.1, 
                               tolerance = 1e-8){

  y <- median(x) #start at median
  converged <- FALSE

  ps <- pseudohuber_loss(x -y, d)
  iter <- 1
  #successive updates using quadratic approx to pseudohuber criterion
  #detailed in supplement of Clausen & Willis (2024)
  while(!converged){
    scaled_sq <- ((x - y)/d)^2

    w <- sqrt(1/(1 + scaled_sq))
    old_y <- y
    y <- weighted.mean(x,w)
    ps <- c(ps,pseudohuber_loss(x- y,d))

    iter <- iter + 1

    if(abs(y - old_y)<tolerance){
      converged <- TRUE
    }

  }

  return(y)
}

#' @rdname pseudohuber_median
#' @export
psuedohuber_median <- pseudohuber_median


