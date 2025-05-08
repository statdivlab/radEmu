
#pseudohuber criterion
pseudohuber_loss <- function(x,#values to computed criterion from 
                        d #smoothing parameter (how to compromise between mean and median)
                        ){
  return(sum(d^2*(sqrt(1 + (x/d)^2) - d)))
}
