# Sum positive numbers via their logarithms
# Given x and y, this function returns log(exp(x) + exp(y)).
sum_of_logs <- function(x, y = NULL){ #calculates log(sum(exp(x)))

  if(is.null(y)){
    if(length(x) ==1){
      return(x)
    } else{
      if(sum(!is.infinite(x)) == 0){
        if(sum(sign(x) != -1) == 0){
          return(-Inf)
        } else{
          return(Inf)
        }
      }
      x <- x[order(x, decreasing = TRUE)]
      return(x[1] + log(1 + sum(exp(x[2:length(x)] - x[1]))))}
  } else{
    if(length(x) != length(y)){
      stop("If y is provided, x and y must be the same length")
    }
    n <- length(x)
    maxes <- pmax(x,y)
    mins <- pmin(x,y)
    return(maxes + log(1 + exp(mins - maxes)))
  }
}
