
micro_wald <- function(Y,
                             X,
                             B,
                             test_kj,
                             constraint_fn,
                             constraint_grad_fn,
                             nominal_coverage = 0.95){
  n <- nrow(Y)
  J <- ncol(Y)
  p <- ncol(X)

  z <- update_z_no_wts(Y,X,B)

  X_cup = X_cup_from_X(X,J)
  scores <- vector(n,mode = "list")

  #compute score contributions of observations i = 1 through n
  for(i in 1:n){
    X_cup_i <- X_cup[(i - 1)*J + 1:J,]
    scores[[i]] <- as.matrix(dpll_dB_cup(X[i,,drop = FALSE],Y[i,,drop = FALSE],B))
  }

  Dy <- Reduce("+",lapply(scores,function(x) tcrossprod(x)))
  B_cup <- B_cup_from_B(B)
  I <- f_info(Y,B_cup,B,X,X_cup)
  I_geninv <- eigen(I)
  n_eigen <- p*(J - 1)
  n_remaining <- p*J - n_eigen
  I_geninv <- I_geninv$vectors%*%diag(c(1/I_geninv$values[1:n_eigen],rep(0,n_remaining)))%*%t(I_geninv$vectors)

  rob_cov <- I_geninv%*%Dy%*%I_geninv

  for(s in 1:nrow(test_kj)){

    null_k <- test_kj$k[s]
    null_j <- test_kj$j[s]

  H <- matrix(0,nrow = p, ncol = J)
  for(j in 1:J){
    H[null_k,j] <- constraint_grad_fn(B[null_k,])[j]
  }
  H[null_k,null_j] <- H[null_k,null_j] -1
  H_cup <- B_cup_from_B(H)

  var_kj <- as.numeric(as.matrix(Matrix::t(H_cup)%*%rob_cov%*%H_cup))*(n/(n - p))
  # print(var_kj)
  test_kj[s,"se"] <- sqrt(var_kj)

  ci <- B[null_k,null_j] + c(-1,1)*qnorm(1 - (1- nominal_coverage)/2)*sqrt(var_kj)
  z <- B[null_k,null_j]/sqrt(var_kj)
  pval <- pchisq(z^2,1,lower.tail = FALSE)
  test_kj[s,c("lower","upper")] <- ci
  test_kj[s,"pval"] <- pval
  }

  return(test_kj)

}


