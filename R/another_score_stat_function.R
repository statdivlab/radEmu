

another_score_stat_function <-
  function(Y,
           X,
           B,
           j_constr,
           k_constr,
           j_ref,
           constraint_fn,
           model_based = FALSE){


    J <- ncol(Y)
    p <- ncol(X)
    n <- nrow(Y)

    for(k in 1:p){
      B[k,] <- B[k,] - B[k,j_ref]
    }

    z <- update_z(Y = Y,
                  X = X,
                  B = B)

    dB <- vector(n, mode = "list")
    ddB <- vector(n,mode = "list")
    log_means <- 0*Y

    ll_fn <- function(B_cup,i){
      B <- B_from_B_cup(B_cup,J,p)
      z <- update_z(Y,X,B)
      logmeans <- sapply(1:J,
                      function(j)
                        X[i,,drop = FALSE]%*%B[,j,drop = FALSE] +
                        z[i])
      return(sum(Y[i,]*logmeans - exp(logmeans)))
    }


    for(i in 1:n){
      dB[[i]] <- numDeriv::grad(function(x) ll_fn(matrix(x,ncol = 1),i),
                                as.numeric(as.matrix(B_cup_from_B(B))))

      ddB[[i]] <- numDeriv::hessian(function(x) ll_fn(matrix(x,ncol = 1),i),
                                 as.numeric(as.matrix(B_cup_from_B(B))))


    }

    I <- -1*Reduce("+",ddB)


    Dy <- Reduce("+",lapply(dB, function(x) x%*%t(x)))
    score <- Reduce("+",dB)
    score <- matrix(score,ncol = 1)

    H <- 0*B
    H[k_constr,] <- numDeriv::grad(function(x) constraint_fn(x) - x[j_constr],B[k_constr,])

    H_cup <- B_cup_from_B(H)

    remove_mat <- 0*B
    remove_mat[,j_ref] <- 1
    remove <- B_cup_from_B(remove_mat)
    remove <- which(remove==1)

    score <- score[-remove,]
    I <- I[-remove,-remove]
    Dy <- Dy[-remove,-remove]
    H_cup <- H_cup[-remove]

    I_inv <- solve(I)
    if(model_based){
      return(
        t(score)%*%I_inv%*%H_cup%*%
          solve(t(H_cup)%*%I_inv%*%H_cup)%*%
          t(H_cup)%*%I_inv%*%score
      )}
      if(!model_based){
        return(
          t(score)%*%I_inv%*%H_cup%*%
            solve(t(H_cup)%*%I_inv%*%Dy%*%I_inv%*%H_cup)%*%
            t(H_cup)%*%I_inv%*%score
        )
      }



  }
