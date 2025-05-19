
#function that does Wald testing / computes Wald CIs
micro_wald <- function(Y,
                       X,
                       X_cup,
                       B,
                       test_kj,
                       constraint_fn,
                       constraint_grad_fn,
                       nominal_coverage = 0.95,
                       verbose = FALSE,
                       return_info = TRUE,
                       return_Dy = TRUE,
                       return_sandwich = FALSE,
                       j_ref = NULL,
                       cluster = NULL){
  n <- nrow(Y)
  J <- ncol(Y)
  p <- ncol(X)
  
  
  if(is.null(j_ref)){
    j_ref <- get_j_ref(Y)
  }
  #impose convenience constraint and update z
  for(k in 1:p){
    B[k,] <- B[k,] - B[k,j_ref]
  }
  z <- update_z(Y,X,B)
  #long form B
  B_cup <- B_cup_from_B(B)
  
  #drop columns corresp. to j_ref from X_cup
  to_erase <- (j_ref - 1)*p + 1:p
  # X_cup_smaller <- X_cup[,-to_erase]
  
  
  
  scores <- vector(n,mode = "list")
  
  
  if(verbose){
    message("Computing 'meat' matrix.")
  }
  
  
  log_means <- X%*%B + matrix(z,ncol = 1)%*%matrix(1,nrow = 1, ncol = J)
  
  #compute residuals
  Y_diff <- Y - exp(log_means)
  
  #compute score from residuals
  scores <- lapply(1:n,
                   function(i){
                     # print(i);
                     Y_diff[i,,drop = FALSE]%*%
                       X_cup[(i - 1)*J + 1:J,]})
  
  if(!is.null(cluster)){
    scores <- lapply(unique(cluster),
                     function(i) 
                       Reduce("+",scores[cluster == i]))
  }

  
  score_mat <- do.call(rbind,scores)
  score_mat <- methods::as(score_mat,"sparseMatrix")
  #empirical score covariance
  Dy <- Matrix::crossprod(score_mat)
  
  
  
  if(verbose){
    message("Computing information matrix.")
  }
  
  #get information matrix
  I <- f_info(Y,B_cup = B_cup_from_B(B),B,X,X_cup,compute_together = FALSE)
  
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
  for(s in 1:nrow(test_kj)){
    
    
    null_k <- test_kj$k[s]
    null_j <- test_kj$j[s]
    
    if(verbose){
      message("Performing test ", s, " of ", nrow(test_kj),": row ", null_k,
              " and column ", null_j, " of B.")
    }
    
    H <- matrix(0,nrow = p, ncol = J - 1 )
    H[null_k,] <- constraint_grad_fn[[null_k]](B[null_k,])[-j_ref]
    
    
    if(null_j != j_ref){
      null_j_index <- ifelse(null_j< j_ref,null_j,null_j - 1)
      H[null_k,null_j_index] <-  H[null_k,null_j_index] - 1
    }
    H_cup <- B_cup_from_B(H)
    
    var_kj <- sum(as.numeric(as.matrix(half_rob_cov%*%H_cup)^2))
    
    test_kj[s,"se"] <- sqrt(var_kj)
    
    ci <- B[null_k,null_j] + c(-1,1)*qnorm(1 - (1- nominal_coverage)/2)*sqrt(var_kj)
    z_stat <- B[null_k,null_j]/sqrt(var_kj)
    pval <- pchisq(z_stat^2,1,lower.tail = FALSE)
    test_kj[s,c("lower","upper")] <- ci
    test_kj[s,"pval"] <- pval
  }
  
  return(list(coefficients = test_kj,
              I = I,
              Dy = Dy))
  
  
}


