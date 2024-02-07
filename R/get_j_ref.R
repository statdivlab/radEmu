get_j_ref <- function(Y){
  return(which.max(colSums(Y>0))) ### which.max returns first index if there are multiple
}