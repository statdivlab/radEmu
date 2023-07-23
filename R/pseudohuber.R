

pseudohuber <- function(x,d){
  return(sum(d^2*(sqrt(1 + (x/d)^2) - d)))
}
