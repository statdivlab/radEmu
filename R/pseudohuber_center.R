
#get the pseudohuber smoothed median
pseudohuber_center <- function(x,
                               d = 1, #smoothing parameter
                               limit = 20, #max iterations
                               tolerance = 1e-8){

  y <- median(x) #start at median
  converged <- FALSE

  ps <- pseudohuber(x -y, d)
  iter <- 1
  #successive updates using quadratic approx to pseudohuber criterion
  #detailed in supplement of Clausen & Willis (2024)
  while(!converged){
    scaled_sq <- ((x - y)/d)^2

    w <- sqrt(1/(1 + scaled_sq))
    old_y <- y
    y <- weighted.mean(x,w)
    ps <- c(ps,pseudohuber(x- y,d))

    iter <- iter + 1

    if(abs(y - old_y)<tolerance){
      converged <- TRUE
    }

  }

  return(y)
}




