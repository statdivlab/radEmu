get_j_ref <- function(Y){
  return(which.max(colSums(Y>0)))
}