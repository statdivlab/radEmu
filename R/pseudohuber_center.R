

pseudohuber_center <- function(x,
                               d = 1,
                               limit = 20,
                               tolerance = 1e-8){

  # return(optim(median(x),
  #       function(y) pseudohuber(x - y,d),
  #       lower = -limit,
  #       upper = limit,
  #       method = "Brent")$par)

  y <- median(x)
  converged <- FALSE

  ps <- pseudohuber(x -y, d)
  iter <- 1
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
#
#     print(y)
#     print(abs(y - old_y))

  }

  return(y)
}




