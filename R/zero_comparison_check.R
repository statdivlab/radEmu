# function to check for zero comparison parameters 
zero_comparison_check <- function(X, Y) {
  
  # remove intercept 
  base_X <- X[, -1, drop = FALSE]
  
  # check if there are more than 1 columns left 
  if (ncol(base_X) > 1) {
    
    # check if there is a single categorical covariate 
    # i.e. all rows add up to 0 or 1 
    if (all(rowSums(base_X) %in% c(0, 1))) {
      
      # get indices for each group 
      n_groups <- ncol(base_X) + 1
      group_ind <- vector(mode = "list", length = n_groups)
      group_ind[[1]] <- which(rowSums(base_X) == 0)
      
      for (i in 1:ncol(base_X)) {
        group_ind[[i + 1]] <- which(base_X[, i] == 1)
      }
      
      # get Y sums for each group 
      J <- ncol(Y)
      group_counts <- matrix(NA, nrow = n_groups, ncol = J)
      for (i in 1:n_groups) {
        group_counts[i, ] <- colSums(Y[group_ind[[i]], ])
      }
      
      # get matrix that is (p - 1) x J that gives whether or not parameter is zero-comparison
      zero_comp <- matrix(NA, nrow = n_groups - 1, J)
      for (i in 1:(n_groups - 1)) {
        zero_comp[i, ] <- (group_counts[1, ] == 0) * (group_counts[i + 1, ] == 0) == 1
      }
      
      # if there are any zero-comparison parameters
      if (sum(zero_comp) > 0) {
        cov <- colnames(base_X)
        cat <- colnames(Y)
        zero_comp_dat <- data.frame(covariate = rep(cov, each = length(cat)),
                                    category = rep(cat, length(cov)))
        zero_comp_dat$zero_comparison <- as.vector(t(zero_comp))
        return(zero_comp_dat)
      }
      
    }
    
  }
  
  # if there are no zero-comparison parameters, return NULL
  return(NULL)
  
}