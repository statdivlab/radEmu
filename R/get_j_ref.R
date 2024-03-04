get_j_ref <- function(Y){
  return(which.max(colSums(Y>0))) ## which.max gives first location in event of ties
}