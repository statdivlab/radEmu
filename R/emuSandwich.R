
emuSandwich <-  function(emuMod,
                         nsim,
                         alpha_level = 0.05){

    rect_weights <- emuMod$weights
    weights <- as.numeric(Y_to_Y_tilde(emuMod$weights))

    X <- emuMod$X
    J <- ncol(emuMod$Y)
    n <- nrow(emuMod$Y)
    p <- ncol(X)
    X_tilde <- X_to_X_tilde(X,J)
    Y <- emuMod$Y
    Y_tilde <- Y_to_Y_tilde(Y)
    S <- Matrix::sparseMatrix(i = 1:(n*J),
                              j = rep(1:n,each = J),
                              x = rep(1, n*J))
    D_tilde <- cbind(X_tilde,S)

    z <- emuMod$z
    B <- emuMod$B

    rownames(D_tilde) <- 1:(n*J)
    B_tilde <- B_to_B_tilde(B)
    theta <- rbind(B_tilde,Matrix::Matrix(z,ncol = 1))
    X_tilde_repar <- X_tilde
    X_tilde_repar <- X_tilde_repar[,1:(p*(J - 1))]
    D_tilde_repar <- cbind(X_tilde_repar,S)
    W <- Matrix::Diagonal(x = as.numeric(exp((D_tilde%*%theta)))/weights)
    W_inv <-  Matrix::Diagonal(x = weights/(
      as.numeric(exp((D_tilde%*%theta)))))
    eig <- eigen(Matrix::crossprod(D_tilde_repar))

    fitted <- exp(X%*%B + matrix(z,ncol = 1)%*%matrix(1, nrow = 1, ncol = J))

    for(j in 1:J){
      B[,j] <- B[,j] - B[,J]
    }

    z <- update_z(Y = Y,
                  rect_weights = rect_weights,
                  X = X,
                  B = B)



    if(emuMod$method=="FL"){
      W_half <- Matrix::Diagonal(x =
                                   sqrt(as.numeric(weights*exp((D_tilde%*%theta)))))
      # info <- Matrix::crossprod(D_tilde_repar,W_half)%*%W_half%*%D_tilde_repar
      info <- Matrix::crossprod(D_tilde_repar,W_half)
      info <- Matrix::tcrossprod(info)

      # message("Computing LU decomposition of information")
      # info_lu <- Matrix::lu(info)
      info_chol <- suppressWarnings(
        try(Matrix::chol(info),silent = TRUE)
      )
      # info_chol <- Matrix::Cholesky(info,pivot = TRUE)
      # info_chol_inv <- Matrix::solve(info_chol)
      if(!inherits(info_chol,"Matrix")){
        perturb <- 1e-8
        while(!inherits(info_chol,"Matrix")){
          # stop("Cholesky decomposition failed")
          message("Reattempting Cholesky decomposition")
          info_chol <- suppressWarnings(
            try(Matrix::chol(info +
                               Matrix::Diagonal(
                                 x = rep(perturb*max(info),
                                         nrow(info)),
                               )),silent = TRUE)
          )
          perturb <- perturb*10
        }
      }
      info_chol_inv <- Matrix::solve(info_chol)


      augmentations <- numeric(nrow(D_tilde))

      for(j in 1:J){
        # message(paste("Hat diagonals for taxon",j,"of",J))
        rel_indices <- sapply(1:n, function(i) (i - 1)*J + j)

        DW <- Matrix::crossprod(D_tilde_repar[rel_indices,],
                                W_half[rel_indices,rel_indices])

        hmm <- Matrix::crossprod(info_chol_inv,DW)


        augmentations[rel_indices] <- #pmax(Matrix::colSums(DW*Z2)/2,0)
          pmax(Matrix::colSums(hmm*hmm)/2,0)

      }
      Y_augmented <- Y_tilde_to_Y(Y_tilde + augmentations,J = J)
      Y <- Y_augmented
    }

    meat <- Matrix::Diagonal(x = rep(1,p*(J-1)+ n))
    bread <- meat

    for(i in 1:n){
      tilde_indices <- (i - 1)*J + 1:J
      D_tilde_i <- D_tilde_repar[tilde_indices,]

      wakeD_i <- Matrix::Diagonal(x = fitted[i,])%*%
                                    D_tilde_i

      W_inv_i <- W_inv[tilde_indices,tilde_indices]

      # var_y_i <-
      #   Matrix::crossprod(Matrix::Matrix( (Y[i,] - fitted[i,]),nrow = 1))
        # Matrix::Diagonal(x =
        #                    (Y[i,] - fitted[i,])^2)

      half_meat <- Matrix::crossprod(wakeD_i,

        W_inv_i%*%Matrix::Matrix( sqrt(.9)*(Y[i,] - fitted[i,]),ncol = 1))

      other_half <-
        Matrix::crossprod(wakeD_i,

                          W_inv_i%*%Matrix::Diagonal(x =  sqrt(.1)*
                                                      sqrt(fitted[i,]/rect_weights[i,])))

      meat <- meat + Matrix::tcrossprod(half_meat) +
        Matrix::tcrossprod(other_half)
      #   Matrix::crossprod(wakeD_i,
      #                     W_inv_i)%*%
      #   var_y_i %*%
      #   W_inv_i%*%wakeD_i
      #
      half_bread <- Matrix::crossprod(wakeD_i,
                                      sqrt(W_inv_i))
      bread <- bread + Matrix::tcrossprod(half_bread)

      # bread <- Matrix::nearPD(bread)$mat
      # meat <- Matrix::nearPD(meat)$mat


    }
    # bread <- bread - Matrix::Diagonal(x = rep(1,p*(J-1)+ n))
    # meat <- meat - Matrix::Diagonal(x = rep(1,p*(J-1)+ n))
    # bread <- Matrix::nearPD(bread)$mat
    # meat <- Matrix::nearPD(meat)$mat
#
#     bread <- Matrix::nearPD(bread)$mat
    # meat <- Matrix::nearPD(meat)$mat
    # #
    # meat_chol <- chol(meat)
    # #
    # cov_half <- qr.solve(bread,meat_chol,tol = 1e-18)



    # bread_inv <- try(solve(bread),silent = TRUE)
    # perturbation <- 1e-12
    # perturb_counter <- 0
    # while(!is.matrix(bread_inv) & perturb_counter < 20){
    # bread_inv <-
    #   try(Matrix::solve(bread +
    #                       Matrix(Diagonal(x =
    #                                         rep(perturbation*max(bread),ncol(bread))))))
    # perturbation <- perturbation*10
    # perturb_counter <- perturb_counter + 1
    # }
    # if(!is.matrix(bread)){
    #   bread_inv <- Matrix::solve(Matrix::nearPD(bread)$mat)
    #   bread_inv <- as.matrix(bread_inv)
    # }
    # if(!is.matrix(bread_inv)){
    #   stop("Inverting bread matrix failed")
    # }
    # covar <- bread_inv%*%meat%*%bread_inv

    eigenmeat <- eigen(meat,symmetric = TRUE)
    half_meat <- eigenmeat$vectors%*%diag(sqrt(pmax(eigenmeat$values,
                                                 0)))
    # eigenbread <- eigen(bread, symmetric = TRUE)
    # half_inv_bread <- eigenbread$vectors%*%diag(sapply(eigenbread$values,
    #                                                    function(x)
    #                                                      ifelse(x>1e-4,
    #                                                             1/sqrt(x),0))) %*%
    #   t(eigenbread$vectors)
    # cov_sqrt <- half_inv_bread%*%half_meat
    # bread_inv <- solve(bread + diag(1e-2*diag(bread)))
    cov_sqrt <- qr.solve(bread ,
                                            half_meat,tol = 1e-20)
    #
    #
    # eigenvar <- base::eigen(covar,symmetric = TRUE)
    #
    # cov_sqrt <- eigenvar$vectors%*%tcrossprod(
    #   diag(sqrt(pmax(0,eigenvar$values))),
    #   eigenvar$vectors)


    constraint_fn <- emuMod$constraint_fn

    sim_draws <- matrix(nrow = J*p,
                        ncol = nsim)
    for(sim in 1:nsim){
      # print(sim)
      presim_b <- rnorm((J - 1)*p)
      presim_z <- rnorm(n)

      cor_sim <- cov_sqrt%*%matrix(c(presim_b,
                                             presim_z),ncol = 1)

      for(k in 1:p){
        b_draw <- c(cor_sim[1:(J - 1) + (k-1)*(J - 1)],0) +
          B[k,]
        b_draw <- b_draw - constraint_fn(b_draw)
        sim_draws[1:J + (k - 1)*J,sim] <- b_draw
      }



    }

    theta <- c(do.call(c,lapply(1:p,
                              function(k) B[k,])),z)

    uppers <- apply(sim_draws,1,function(x) quantile(x,1 - alpha_level/2))
    lowers <- apply(sim_draws,1,function(x) quantile(x,alpha_level/2))
    pvals <- do.call(c,lapply(1:nrow(sim_draws),
    function(k) mean(abs(
      sim_draws[k,] - theta[k]) >= abs(theta[k]))))

    return(list("uppers" = uppers,
                "lowers" = lowers,
                "pvals" = pvals))
}
