
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
                       j_ref = NULL){
  n <- nrow(Y)
  J <- ncol(Y)
  p <- ncol(X)


  if(is.null(j_ref)){
    nice_ref <- which.max(colSums(Y>0))
  } else{
    nice_ref <- j_ref
  }
  for(k in 1:p){
    B[k,] <- B[k,] - B[k,nice_ref]
  }
  z <- update_z(Y,X,B)
  B_cup <- B_cup_from_B(B)

  #drop columns corresp. to nice_ref from X_cup
  to_erase <- (nice_ref - 1)*p + 1:p
  # X_cup_smaller <- X_cup[,-to_erase]



  scores <- vector(n,mode = "list")


  if(verbose){
    message("Computing 'meat' matrix.")
  }

  #compute score contributions of observations i = 1 through n

  if(FALSE){
  for(i in 1:n){
    if(verbose){
      message("Computing score component ", i, " of ", n,".")
    }
    # print(i)
    X_cup_i <- X_cup[(i - 1)*J + 1:J,]
    # scores[[i]] <- as.matrix(dpll_dB_cup(X[i,,drop = FALSE],Y[i,,drop = FALSE],B))
    log_means_i <- X_cup_i%*%B_cup + z[i]
    scores[[i]] <- Matrix::crossprod(X_cup_i,t(Y[i,,drop = FALSE]) - exp(log_means_i))
  }


  } else{
    # score_mat <- methods::as(matrix(0,nrow = n, ncol = ncol(X_cup)),"sparseMatrix")
    # Matrix::sparseMatrix()
    #
    log_means <- X%*%B + matrix(z,ncol = 1)%*%matrix(1,nrow = 1, ncol = J)
    #
    # for(i in 1:n){
    #   if(verbose){
    #     message("Computing score component ", i, " of ", n,".")
    #   }
    #   # print(i)
    #   # X_cup_i <- X_cup[(i - 1)*J + 1:J,]
    #   # # scores[[i]] <- as.matrix(dpll_dB_cup(X[i,,drop = FALSE],Y[i,,drop = FALSE],B))
    #   # log_means_i <- X_cup_i%*%B_cup + z[i]
    #   X_i <- X[i,]
    #   diffs <- Y[i,] - exp(log_means[i,])
    #   score_mat[i,] <- do.call(c,
    #                            lapply(1:J,
    #                            function(j)
    #                            diffs[j]*(X_i)))
    # }

    Y_diff <- Y - exp(log_means)

    scores <- lapply(1:n,
                     function(i){
                       # print(i);
                       Y_diff[i,,drop = FALSE]%*%
                                                   X_cup[(i - 1)*J + 1:J,]})
    # Dy <- Matrix::crossprod(scores[[1]])
    # for(i in 2:n){
    #   # if(verbose){
    #   #   message("Computing meat component ", i, " of ", n,".")
    #   # }
    #   print(i)
    #   Dy <- Dy + Matrix::crossprod(scores[[i]])
    # }
    score_mat <- do.call(rbind,scores)
    score_mat <- methods::as(score_mat,"sparseMatrix")
    Dy <- Matrix::crossprod(score_mat)


    # message("Computing 'meat' matrix")

    # Dy <- Matrix::crossprod(score_mat)

  }

  if(verbose){
    message("Computing information matrix.")
  }

  # if(ncol(X_cup)<2500){
  I <- f_info(Y,B_cup = B_cup_from_B(B),B,X,X_cup,compute_together = FALSE)

  if(verbose){
    message("Inverting information to obtain 'bread' matrix.")
  }

  half_rob_cov <- Matrix::t(Matrix::solve(I[-to_erase,-to_erase], Matrix::t(score_mat)[-to_erase,],
                                method = "cholmod_solve"))

  #return to original parametrization
  for(k in 1:p){
    B[k,] <- B[k,] - constraint_fn(B[k,])
  }

  if(verbose){
    message("Performing Wald tests and constructing Wald CIs.")
  }
  for(s in 1:nrow(test_kj)){



    null_k <- test_kj$k[s]
    null_j <- test_kj$j[s]

    if(verbose){
      message("Performing test ", s, " of ", nrow(test_kj),": row ", null_k,
              " and column ", null_j, " of B.")
    }

  H <- matrix(0,nrow = p, ncol = J - 1 )
  H[null_k,] <- constraint_grad_fn(B[null_k,])[-nice_ref]


  if(null_j != nice_ref){
    null_j_index <- ifelse(null_j< nice_ref,null_j,null_j - 1)
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


