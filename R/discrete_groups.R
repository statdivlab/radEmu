#group rows of a design matrix by covariate pattern
#returns the distinct rows of X (in order of first appearance) and, for each
#row of X, the index of its pattern among the distinct rows
discrete_groups <- function(X){
  key <- apply(X, 1, paste, collapse = "_")
  first <- !duplicated(key)
  list(distinct_X = X[first, , drop = FALSE],
       group = match(key, key[first]))
}

