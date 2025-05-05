#' Test a linear combination
#' 
#' Test a hypothesis about a linear combination of parameters. Currently this is 
#' only implemented for robust Wald tests. 
#'
#' @param fitted_model An \code{emuFit} fitted model. 
#' @param linear_combo A vector of length \code{p = ncol(X)} to determine the linear 
#' combination of parameters to be tested. 
#' @param hypothesis The null hypothesis to compare the linear combination of
#' coefficients against. The default value is \code{0}.
#' @param j A vector of \code{j} values giving the indices of taxa to test. The default
#' behavior is to include all taxa. 
#' 
#' @return A data frame giving confidence intervals using robust standard errors, test
#' statistics for the robust Wald test, and p-values, for the linear combination of 
#' parameters that is being tested, for each taxon specified in \code{j}. 
#' 
#' @examples
#' data(wirbel_sample_small)
#' data(wirbel_otu_small)
#' emuRes <- emuFit(formula = ~ Group + Country, data = wirbel_sample_small, Y = wirbel_otu_small,
#'                  test_kj = data.frame(k = 2, j = 1), tolerance = 0.01, 
#'                  constraint_tol = 0.01, B_null_tol = 0.01) 
#'  # here we set large tolerances for the example to run quickly, 
#'  # but we recommend smaller tolerances in practice
#' # test whether CountryFRA - constraint_FRA = CountryUSA - constraint_USA for taxon 1 
#' linear_combo_res <- test_linear_combo(fitted_res = emuRes, 
#'                                      linear_combo = c(0, 0, 1, -1),
#'                                      j = 1)
#'
#' @export
#'
test_linear_combo <- function(fitted_model,
                              linear_combo,
                              hypothesis = 0,
                              j = NULL) {
  
  # get useful arguments
  X <- fitted_model$X
  if (is.null(fitted_model$Y_augmented)) {
    Y <- fitted_model$Y
  } else {
    Y <- fitted_model$Y_augmented
  }
  B <- fitted_model$B
  constraint_fn <- fitted_model$constraint_fn
  constraint_grad_fn <- fitted_model$constraint_grad_fn
  cluster <- fitted_model$cluster
  call_args <- as.list(fitted_model$call)[-1]
  if ("verbose" %in% names(call_args)) {
    verbose <- call_args$verbose %in% c(TRUE, "development")
  } else {
    verbose <- FALSE
  }
  if ("alpha" %in% names(call_args)) {
    alpha <- call_args$alpha
  } else {
    alpha <- 0.05
  }
  nominal_coverage <- 1 - alpha
  
  # set j if needed 
  n <- nrow(X)
  p <- ncol(X)
  J <- ncol(B)
  if (is.null(j)) {
    j <- 1:J
  }
  
  # calculate other needed quantities 
  j_ref <- get_j_ref(Y)
  X_cup <- X_cup_from_X(X, J)
  coefs <- data.frame(j = j, linear_combo_estimate = NA, ci_lower = NA, ci_upper = NA,
                      test_stat = NA, wald_pval = NA)
  
  # copy wald testing from micro_wald, with updates for linear combo 
  
  # impose convenience constraint and update z
  for (k in 1:p) {
    B[k, ] <- B[k, ] - B[k, j_ref]
  }
  z <- update_z(Y, X, B)
  # long form B
  B_cup <- B_cup_from_B(B)
  
  # drop columns corresp. to j_ref from X_cup
  to_erase <- (j_ref - 1) * p + 1:p
  scores <- vector(n, mode = "list")

  if(verbose){
    message("Computing 'meat' matrix.")
  }
  
  log_means <- X%*%B + matrix(z,ncol = 1)%*%matrix(1,nrow = 1, ncol = J)
  
  #compute residuals
  Y_diff <- Y - exp(log_means)
  
  #compute score from residuals
  scores <- lapply(1:n,
                   function(i){
                     Y_diff[i,,drop = FALSE]%*%
                       X_cup[(i - 1)*J + 1:J,]})
  
  if(!is.null(cluster)){
    scores <- lapply(unique(cluster),
                     function(i) 
                       Reduce("+",scores[cluster == i]))
  }
  
  score_mat <- do.call(rbind,scores)
  score_mat <- methods::as(score_mat,"sparseMatrix")
  
  if(verbose){
    message("Computing information matrix.")
  }
  
  # get information matrix
  if (is.null(fitted_model$I)) {
    I <- f_info(Y,B_cup = B_cup_from_B(B), B, X, X_cup, compute_together = FALSE)
  } else {
    I <- fitted_model$I
  }
  
  if(verbose){
    message("Inverting information to obtain 'bread' matrix.")
  }
  
  #get a matrix whose crossproduct with itself is the sandwich covariance for B
  half_rob_cov <- Matrix::t(Matrix::solve(I[-to_erase,-to_erase], Matrix::t(score_mat)[-to_erase,],
                                          method = "cholmod_solve"))
  
  #return to original parametrization
  for(k in 1:p){
    B[k,] <- B[k,] - constraint_fn[[k]](B[k,])
  }
  
  if(verbose){
    message("Performing Wald tests and constructing Wald CIs.")
  }
  
  #testing / confidence interval construction
  for(s in 1:nrow(coefs)){
    
    null_j <- coefs$j[s]
    
    if(verbose){
      message("Performing test ", s, " of ", nrow(coefs),": column ", null_j, " of B.")
    }
    
    H <- matrix(0,nrow = p, ncol = J - 1 )
    
    for (k in 1:p) {
      if (!(linear_combo[k] == 0)) {
        H[k,] <- constraint_grad_fn[[k]](B[k, ])[-j_ref]
        if(null_j != j_ref){
          null_j_index <- ifelse(null_j < j_ref, null_j, null_j - 1)
          H[k, null_j_index] <-  H[k, null_j_index] - linear_combo[k]
        }
      }
    }
    
    H_cup <- B_cup_from_B(H)
    
    var_kj <- sum(as.numeric(as.matrix(half_rob_cov%*%H_cup)^2))
    
    h <- as.vector(B[, null_j] %*% linear_combo - hypothesis)
    coefs$linear_combo_estimate[s] <- h
    ci <- h + c(-1,1)*qnorm(1 - (1- nominal_coverage)/2)*sqrt(var_kj)
    coefs[s, c("ci_lower", "ci_upper")] <- ci
    test_stat <- h/sqrt(var_kj)
    pval <- pchisq(test_stat^2, 1, lower.tail = FALSE)
    coefs[s, c("test_stat", "wald_pval")] <- c(test_stat, pval)
  }
  
  return(coefs)
}