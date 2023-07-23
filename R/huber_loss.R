
huber_loss <- function(x,
                       d = 1){

  loss <- 0*x
  above_cutoff <- abs(x)>d
  loss[above_cutoff] <- d*(abs(x[above_cutoff]) - 0.5*d)
  loss[!above_cutoff] <- 0.5*(x[!above_cutoff])^2
  return(sum(loss))
}
