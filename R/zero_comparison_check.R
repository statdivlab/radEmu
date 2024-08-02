# function to check for zero comparison parameters 
zero_comparison_check <- function(X, Y) {
  
  # X is a matrix from model.matrix and has "assign" attribute 
  if ("assign" %in% names(attributes(X))) {
    
    col_assign <- attr(X, "assign")
    col_shared <- sapply(col_assign, function(x) {sum(col_assign == x)})
    col_kept <- col_assign[col_shared > 1]
    
    if (max(col_shared) > 1) {
      
      zero_comp_dat <- data.frame(covariate = NULL, category = NULL, 
                                  zero_comparison = NULL)
      
      cat_covs <- unique(col_kept)
      for (col in cat_covs) {
        base_X <- X[, col_assign == col]
        
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
          new_zero_comp_dat <- data.frame(covariate = rep(cov, each = length(cat)),
                                      category = rep(cat, length(cov)))
          new_zero_comp_dat$zero_comparison <- as.vector(t(zero_comp))
          zero_comp_dat <- rbind(zero_comp_dat, new_zero_comp_dat)
        }
      }
      if (nrow(zero_comp_dat) > 0) {
        return(zero_comp_dat)
      }
    }
  } else {
  # X is a matrix (not from model matrix), need to check manually
  # this will only identify a singular categorical covariate 
    
    # remove intercept 
    base_X <- X[, -1, drop = FALSE]
    
    # remove any columns of X that include values other than 0 and 1 
    non_cat_cols <- which(apply(base_X, 2, function(x) {sum(!(x %in% c(0, 1)))} > 0))
    if (length(non_cat_cols) > 0 & length(non_cat_cols) == ncol(base_X)) {
      return(NULL)
    } else {
      if (length(non_cat_cols) > 0) {
        base_X <- base_X[, -non_cat_cols, drop = FALSE]
      }
    }
    
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
  }
  
  # if there are no zero-comparison parameters, return NULL
  return(NULL)
  
}