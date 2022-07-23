

#' Fit partially identified log-linear model to 'relative abundance' data
#'
#' This function fits a radEmu model to multivariate outcome data Y. Predictors
#' are specified either via a formula, in which case covariate data must be
#' included as a data frame, or via an already formed model matrix X.
#'
#'
#' @param formula_rhs The right-hand side of a formula specifying which
#' predictors should be included in model
#' @param Y A matrix of outcome data with rows corresponding to samples and
#' columns corresponding to outcome categories (e.g., microbial taxon)
#' @param constraint_fn The function to be used as a constraint -- namely,
#' we enforce constraint_fn(B_k) = 0 for all rows k = 1, ..., p of B. This
#' defaults to median().
#' @param reweight Refit model using weights derived from estimated mean-
#' variance relationship? Default is FALSE.
#
#' @return \item{emuMod}{An emuMod object specifying the model fit and providing
#' point estimates for effects included in model.}
#' @author David Clausen
emuFit <-  function(formula_rhs = NULL,
                    Y,
                    X = NULL,
                    covariate_data = NULL,
                    B = NULL,
                    B_cutoff = 20,
                    tolerance = 1e-1,
                    maxit = 100,
                    verbose = TRUE,
                    constraint_fn = NULL,
                    maxit_glm = 100,
                    optim_only = TRUE,
                    method = "ML",
                    guarded =TRUE,
                    linesearch = FALSE, #only for ML; FL automatically does linesearch
                    reweight = FALSE,
                    reweight_blocks = NULL,
                    weights = NULL,
                    test_firth = FALSE,
                    return_a_lot = TRUE,
                    prefit = TRUE){

  if(!is.null(formula_rhs)){
    if(is.null(covariate_data)){
      stop("If formula_rhs is provided, covariates named in formula must be
           provided inside covariate_data.")
    }
    if(!is.data.frame(covariate_data)){
      stop("Argument covariate_data must be a data frame.")
    }
    if(!inherits(formula_rhs,"formula")){
      formula_rhs <- as.formula(formula_rhs)
    }
    X <- model.matrix(formula_rhs,covariate_data)
  }



  if(is.null(constraint_fn)){
    constraint_fn <- function(x){ median(x)}
  }

  p <- ncol(X)
  if(p<2){
    stop("X must contain an intercept column and at least one other (linearly
independent) column")
  }
  J <- ncol(Y)
  n <- nrow(X)
  if(is.null(B)){
    B <- matrix(0, nrow = p,ncol = J)
  }

  if(reweight){
    if(is.null(B) | prefit){
      if(verbose){
        message("Fitting initial model")
      }
      if(!prefit){
        warning("Ignoring prefit = FALSE; provide value of B to be
used to estimate weights if you wish to skip the initial fitting step.")
      }
      prefit <- emuFit(X = X,
                       Y = Y,
                       B = B,
                       B_cutoff = B_cutoff,
                       tolerance = tolerance,
                       maxit = maxit,
                       constraint_fn = constraint_fn,
                       verbose = verbose,
                       method = method,
                       reweight = FALSE,
                       weights = weights
      )
    } else{
      message("Computing fitted values for reweighting
using input B; to estimate B prior to reweighting, set prefit = TRUE.")
      if(is.null(weights)){
        weights <- 0*Y + 1
      }
      prefit <- list("B" = B,
                     "z" = update_z(Y, weights, X, B))
    }

    prefitted <- exp(X%*%prefit$B +
                       matrix(prefit$z,ncol = 1)%*%matrix(1,nrow = 1, ncol = J))

    if(is.null(reweight_blocks)){
    weights <- prefitted_to_weights(Y = Y,
                                    prefitted = prefitted)
    } else{
      weights <- 1 + 0*Y
      rw_block_uni <- unique(reweight_blocks)

      for(rw_block in rw_block_uni){
        which_rw <- which(reweight_blocks == rw_block)
        weights[which_rw,] <-
          prefitted_to_weights(Y = Y[which_rw,],
                               prefitted = prefitted[which_rw,])
      }
    }


  }

  if(is.null(weights)){
    weights <- rep(1,n*J)
    rect_weights <- 1.0 + 0.0*Y
  } else{
    rect_weights <- weights
    weights <- as.numeric(Y_to_Y_tilde(weights))
  }

  if(method == "FL"){

    X_tilde <- X_to_X_tilde(X,J)
    Y_tilde <- Y_to_Y_tilde(Y)
    S <- Matrix::sparseMatrix(i = 1:(n*J),
                              j = rep(1:n,each = J),
                              x = rep(1, n*J))
    D_tilde <- cbind(X_tilde,S)

    z <- apply(Y*rect_weights,1,function(x) log(sum(x))) -
      apply(exp(X%*%B)*rect_weights,1,function(x) log(sum(x)))

    rownames(D_tilde) <- 1:(n*J)
    B_tilde <- B_to_B_tilde(B)
    theta <- rbind(B_tilde,Matrix::Matrix(z,ncol = 1))
    X_tilde_repar <- X_tilde
    X_tilde_repar <- X_tilde_repar[,-((2:p)*J)]
    D_tilde_repar <- cbind(X_tilde_repar,S[,-n])

    B_tilde <- B_to_B_tilde(B)
    theta <- rbind(B_tilde,Matrix::Matrix(z,ncol = 1))
    W <- Matrix::Diagonal(x = weights*as.numeric(exp((D_tilde%*%theta))))

    lls <- log_likelihood_wide(Y,
                               rect_weights,
                               X,
                               B,
                               z) +
      calculate_firth_penalty(D_tilde = D_tilde,
                              W = W,
                              n_skip = p)


    converged <- FALSE
    lilstep <- 0
    iter <- 1
    while(!converged){
      if(verbose){
        message(paste("Iteration ", iter,"; penalized log likelihood ", lls[iter],
                      sep ="", collapse = ""))
      }

      B_tilde <- B_to_B_tilde(B)
      theta <- rbind(B_tilde,Matrix::Matrix(z,ncol = 1))
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
      # WD_qr <- Matrix::qr(W_half%*%D_tilde)
      # info_lu <- Matrix::expand(info_lu)
      # wx_svd <- svd(W_half%*%D_tilde_repar)

      # tergit <- crossprod(crossprod(
      #   D_tilde_repar,
      #   W_half
      # ),solve(info))%*%
      #                crossprod(
      #                  D_tilde_repar,
      #                  W_half
      #                )
      #
      #
      #
      # terg_aug <- diag(as.matrix(tergit))


      # message("Inverting L and U")
      #
      #     Z1 <- Matrix::solve(a = info_lu$P%*%info_lu$L,
      #                         b = Matrix::t(D_tilde_repar)%*%W_half,
      #                         sparse = FALSE)
      #
      #     Z2 <- Matrix::solve(
      #       a = Matrix::t(info_lu$U%*%info_lu$Q),
      #       b = Matrix::t(D_tilde_repar)%*%W_half,
      #       sparse = FALSE)
      #
      #     augmentations <- Matrix::colSums(Z1*Z2)/2



      # info_eigen <- eigen(info, symmetric = TRUE)
      # info_V <- Matrix::Matrix(info_eigen$vectors)
      # info_inv <- info_V%*%Matrix::tcrossprod(
      #   Matrix::Diagonal(x = sapply(info_eigen$values,
      #                               function(x) ifelse(x>1e-6,1/x,0))),
      #   info_V
      # )
      #
      # info_half_inv <- info_V%*%Matrix::tcrossprod(
      #   Matrix::Diagonal(x = sapply(info_eigen$values,
      #                               function(x) ifelse(x>1e-6,1/sqrt(x),0))),
      #   info_V
      # )

      augmentations <- numeric(nrow(D_tilde))
      # L_inv <- Matrix::solve(Matrix::t(info_lu$P)%*%info_lu$L)
      # U_inv <- Matrix::solve(info_lu$U%*%info_lu$Q)

      # max(abs(with(info_lu,t(P)%*%L%*%U%*%Q) - info))

      # info_inv <- with(info_lu,
      #                  Matrix::solve(info_lu@U[,info_lu@q + 1])%*%
      #                    Matrix::solve(info_lu@L[info_lu@p + 1,])
      # )
      if(verbose){ message("Computing data augmentation. This may take a moment on
larger datasets.")}
      for(j in 1:J){
        # message(paste("Hat diagonals for taxon",j,"of",J))
        rel_indices <- sapply(1:n, function(i) (i - 1)*J + j)

        # Z1 <- Matrix::solve(a = info_lu$P%*%info_lu$L,
        #                     b = Matrix::t(D_tilde_repar[rel_indices,])%*%
        #                       W_half[rel_indices,
        #                              rel_indices],
        #                     sparse = FALSE)
        DW <- Matrix::crossprod(D_tilde_repar[rel_indices,],
                                W_half[rel_indices,rel_indices])
        # Z1 <- L_inv%*%DW
        # Z2 <- U_inv%*%DW
        # Z2 <- info_inv%*%DW
        hmm <- Matrix::crossprod(info_chol_inv,DW)

        # info_inv_alt <- get_info_inv(X_tilde_repar,S,J,n,weights)



        # Z2 <- Matrix::solve(
        #   a = Matrix::t(info_lu$U%*%info_lu$Q),
        #   b = Matrix::t(D_tilde_repar[rel_indices,])%*%
        #     W_half[rel_indices,rel_indices],
        #   sparse = FALSE)


        augmentations[rel_indices] <- #pmax(Matrix::colSums(DW*Z2)/2,0)
          pmax(Matrix::colSums(hmm*hmm)/2,0)
        #   # augmentations[rel_indices] <-
        #   # Matrix::diag(
        #   #   W_half[rel_indices,
        #   #          rel_indices]%*%D_tilde_repar[rel_indices,,drop = FALSE] %*%
        #   #     Matrix::tcrossprod(info_inv + Matrix::Diagonal(x = rep(1e-5,nrow(info_inv))),
        #   #                      W_half[rel_indices,rel_indices]%*%
        #   #                        D_tilde_repar[rel_indices,,drop = FALSE]))
        #   starter <-  Matrix::crossprod(D_tilde_repar[rel_indices,,drop = FALSE],
        #                                 W_half[rel_indices,rel_indices])
        #   seconds <- Matrix::solve(info,
        #                   starter)
        #     #
        #     # W_half[rel_indices,
        #     #                       rel_indices]%*%D_tilde_repar[rel_indices,,drop = FALSE] %*%
        #     # info_half_inv
        #
        #   augmentations[rel_indices] <- colSums(as.matrix(starter)*
        #                                           (as.matrix(seconds)))/2
        #
      }

      Y_augmented <- Y_tilde_to_Y(Y_tilde + augmentations,J = J)
      curr_der <-
        emuDeriv(D_tilde = D_tilde,
                 theta = theta,
                 Y_tilde = Y_tilde + augmentations,
                 weights = weights)

      # ll_fn <- function(theta){
      #   theta <- Matrix(theta,ncol= 1)
      #   sum(Matrix::Diagonal(x = weights)%*%((Y_tilde + augmentations)*(D_tilde%*%theta) -
      #                                      exp(D_tilde%*%theta)))
      # }
      #
      #
      # prop_ll <- log_likelihood_wide(Y = Y_augmented,
      #                                rect_weights = rect_weights,
      #                                X = X,
      #                                B = B_tilde_to_B(theta[,1:20,
      #                                                       drop = FALSE],J = 10),
      #                                z = prop_z) +
      #   calculate_firth_penalty(D_tilde_repar = D_tilde,
      #                           W = prop_W,
      #                           n_skip = p)
      # num_der <- numDeriv::grad(ll_fn, as.numeric(theta))

      message(paste("Norm of gradient",
                    round(sqrt(sum(as.numeric(as.matrix(curr_der))^2)),2)))

      # plot(asinh(curr_der)[order(as.numeric(as.matrix(curr_der)))],pch = ".")
      # points(as.numeric(matrix(theta))[order(as.numeric(as.matrix(curr_der)))],pch = ".",col = "red")
      # message(which.max(abs(as.numeric(as.matrix(curr_der)))))

      if(test_firth){
        augmentations_also <- diag(
          as.matrix(W_half%*%D_tilde_repar%*%Matrix::solve(info)%*%
                      Matrix::t(W_half%*%D_tilde_repar))
        )/2
        return(list("augmentations_chol" = augmentations,
                    "augmentations_naive" = augmentations_also))
      }
      if(verbose){
      message("Computing ML fit on augmented data")}
      ml_fit <- emuFit(X = X,
                       Y = Y_augmented,
                       B = B,
                       B_cutoff = B_cutoff,
                       tolerance = tolerance,
                       maxit = 1,
                       verbose = FALSE,
                       # verbose = TRUE,
                       maxit_glm = maxit_glm,
                       constraint_fn = constraint_fn,
                       method = "ML",
                       reweight = FALSE,
                       weights = rect_weights)

      # message(paste("Norm of update step",
      #               round(sqrt(sum((ml_fit$B - B)^2)),4)))

      B_direction <- ml_fit$B - B
      z_direction <- ml_fit$z - z
      stepsize <- 2

      accepted <- FALSE
      while(!accepted){
        stepsize <- stepsize/2
        prop_B <- B + stepsize*B_direction
        for(k in 1:p){
          prop_B[k,] <- prop_B[k,] - constraint_fn(prop_B[k,])
        }
        prop_z <- z + stepsize*z_direction
      #
      ##################
      # prop_B <- ml_fit$B
      # prop_z <- ml_fit$z
      ##################
        prop_B_tilde <- B_to_B_tilde(prop_B)
        prop_theta <- rbind(prop_B_tilde,Matrix::Matrix(prop_z,ncol = 1))
        prop_W <- Matrix::Diagonal(x = weights*as.numeric(exp((D_tilde%*%prop_theta))))
      #

        firth_pen <-  calculate_firth_penalty(D_tilde = D_tilde,
                                              W = prop_W,
                                              n_skip = p)
        prop_ll <- log_likelihood_wide(Y = Y,
                                       rect_weights = rect_weights,
                                       X = X,
                                       B = prop_B,
                                       z = prop_z) + firth_pen

        # message(paste("Firth penalty term",firth_pen))
      #
        if((prop_ll >= lls[iter])|stepsize <= 1e-6){
          accepted <- TRUE
        }
        accepted <- TRUE
      }
      #
      lls <- c(lls,
               prop_ll)
      B <- prop_B
      z <- prop_z

      # plot(asinh(as.numeric(B)),cex = .3)

      # plot(do.call(c,lapply(1:nrow(B), function(k) B[k,])))

      stepsize <- 1
      if(stepsize ==1){
        lilstep <- 0
      } else{
        lilstep <- lilstep + 1
      }

      step_criterion <- ifelse(stepsize ==1,TRUE,
                               lilstep >2)

      if((abs(lls[iter + 1] - lls[iter]) < tolerance)#&
         # (max(lls[1:iter] - lls[iter + 1]) < tolerance)
          &
         step_criterion
      ){
        converged <- TRUE
      }

      if(iter>=maxit){
        converged <- TRUE
      }
      iter <- iter + 1
    }

    if(return_a_lot){
      return(list("B" = B,
                  "z" = z,
                  "ll" = lls[iter - 1],
                  "Y" = Y,
                  "X" = X,
                  "tolerance" = tolerance,
                  "maxit" = maxit,
                  "constraint_fn" = constraint_fn,
                  "maxit_glm" = maxit_glm,
                  "method" = method,
                  "weights" = rect_weights,
                  "reweight" = reweight))
    } else{
      return(
        list("B" = B,
             "z" = z,
             "ll" = lls[iter - 1],
             "method" = method,
             "weights" = rect_weights,
             "reweight" = reweight)
      )

    }

  }

  z <- apply(Y*rect_weights,1,function(x) log(sum(x))) -
    apply(exp(X%*%B)*rect_weights,1,function(x) log(sum(x)))

  lls <- numeric(maxit + 1)

  lls[1] <- log_likelihood_wide(Y = Y,
                                rect_weights = rect_weights,
                                X = X,
                                B = B,
                                z = z)

  converged <- FALSE
  iter <- 1
  while(!converged){
    if(verbose){
      message(paste("Iteration ", iter,"; log likelihood ", lls[iter],
                    sep ="", collapse = ""))
    }

    B_old <- B
    z_old <- z


    # update_order <- sample(1:J)
    update_order <- 1:J
    for(j in update_order){
      # message(j)

      if(!optim_only){
      Bj_update <- emuFit_one(Y = Y,
                              X = X,
                              rect_weights = rect_weights,
                              j = j,
                              B = B,
                              z = z,
                              method = method,
                              maxit_glm = maxit_glm,
                              WD = NULL,
                              info_inv = NULL)

      if(is.null(Bj_update)){
        stop("glm update failed!")
      }
      #
      # Bj_update[abs(Bj_update)>B_cutoff] <-
      #   B_cutoff*sign(Bj_update[abs(Bj_update)>B_cutoff])
      #
      # step_dir <-  Bj_update - B[,j]
      old_lj <- sum(
      (Y[,j]*(X%*%B[,j] + z) - exp(X%*%B[,j] + z))*rect_weights[,j])

      new_lj <- sum(
        (Y[,j]*(X%*%Bj_update  + z) - exp(X%*%Bj_update + z))*rect_weights[,j])
      } else{
        new_lj <- 0
        old_lj <- 1
      }
      #
      if(new_lj + 1e-4 < old_lj){
        if(verbose & !optim_only){
        message(paste("Update for taxon", j, "does not increase likelihood;",
                      "recomputing with optim."))
          }

        opt_fun <- function(b){
          -1*sum(
          (Y[,j]*(X%*%matrix(b,ncol = 1)  + z) -
             exp(X%*%matrix(b,ncol = 1) + z))*rect_weights[,j])
          }

        gr_fun <- function(b){
          -1*as.numeric(crossprod(X,rect_weights[,j]*(Y[,j] -
             exp(X%*%matrix(b,ncol = 1) + z))))}


        opt_result <- try(optim(B[,j],
                            opt_fun,
                            gr = gr_fun,
                            method = "BFGS"),
                          silent = TRUE)

        if(inherits(opt_result,"try-error")){
          opt_result <- try(optim(B[,j],
                                  opt_fun,
                                  # gr = gr_fun,
                                  method = "BFGS"))
        }

        B[,j] <- opt_result$par
        # message(paste(paste("Beta", j, "updated to"), paste(round(B[,j],2),
        #                                                     collapse = " "),
        #               collapse = ""))

      } else{
        B[,j] <- Bj_update
      }

      # new_lj <- old_lj - 1e4


      # if(!is.null(Bj_update)){
      # linesearch_counter <- 1
      # stepsize <- 1
      # while((new_lj < old_lj)& linesearch_counter<8){
      #   prop_Bj <- B[,j] + step_dir*stepsize
        # new_lj <- sum(
        #   (Y[,j]*(X%*%prop_Bj  + z) - exp(X%*%prop_Bj + z))*rect_weights[,j])
#
#
#
#       if(new_lj >= old_lj){
#       B[,j] <- prop_Bj} else{
#         # stop()
#         # if(j == 760){
#         message(paste("Step too large for outcome ",
#                       j, "; setting step size to ",
#                       stepsize/2,
#                       collapse = "",
#                       sep = ""))
#         # }
#       }
#         linesearch_counter <- linesearch_counter + 1
#         stepsize <- stepsize/2
#       }
#
#       if(new_lj + 1e-4 < old_lj){
#         # stop()s
#       # } else{
#         # message(paste("Skipping update for outcome ", j, "; glm fit failed.",
#         #               sep = "",collapse = ""))
#       }

      z <- update_z(Y = Y,
                    rect_weights = rect_weights,
                    X = X,B)

    }

    for(k in 1:p){
      B[k,] <- B[k,] - constraint_fn(B[k,])
    }

    z <- log(Matrix::rowSums(Y*rect_weights)) -
      log(Matrix::rowSums(exp(X%*%B)*rect_weights))

    stepsize <- 1
    # if(linesearch){
    # B_direction <- B - B_old
    # z_direction <- z - z_old
    # stepsize <- 2
    #
    # accepted <- FALSE
    # while(!accepted){
    #   stepsize <- stepsize/2
    #   message(stepsize)
    #   prop_B <- B_old + stepsize*B_direction
    #
    #   for(k in 1:p){
    #     prop_B[k,] <- prop_B[k,] - constraint_fn(prop_B[k,])
    #   }
    #
    #   prop_z <- log(Matrix::rowSums(Y*rect_weights)) -
    #     log(Matrix::rowSums(exp(X%*%prop_B)*rect_weights))
    #
      # prop_ll <- log_likelihood_wide(Y = Y,
      #                                rect_weights = rect_weights,
      #                                X = X,
      #                                B = prop_B,
      #                                z = prop_z)

      # B_tilde_current <- B_to_B_tilde(B)
      #
      # ll_fn <-function(x){ log_likelihood_wide(Y = Y,
      #                                               rect_weights = rect_weights,
      #                                               X = X,
      #                                               B = B_tilde_to_B(x,
      #                                                                p = p,
      #                                                                J = J),
      #                                               z = z)}
      # #
      # nd <- numDeriv::grad(ll_fn,
      #                      x = as.numeric(B_to_B_tilde(B)))
      #
      # X_tilde <- X_to_X_tilde(X,J)
      # S <- Matrix::sparseMatrix(i = 1:(n*J),
      #                           j = rep(1:n,each = J),
      #                           x = rep(1, n*J))
      # D_tilde <- cbind(X_tilde,S)
      # B_tilde <- B_to_B_tilde(B)
      # theta <- rbind(B_tilde,matrix(z,ncol = 1))
      # Y_tilde <- Y_to_Y_tilde(Y)
      # plot(log(Y_tilde), D_tilde%*%theta)
      #
      # plot(nd, as.matrix(crossprod(Y_tilde - exp(D_tilde%*%theta),D_tilde))[,1:(p*J)])
      # abline(a = 0, b= 1)
      #
      # nd_wide <-
    #
    #   if(prop_ll >= lls[iter]){
    #     accepted <- TRUE
    #   }
    # }
    #
    # lls <- c(lls,
    #          prop_ll)
    #
    # B <- prop_B
    # z <- prop_z
    # }




    # log_means <- X%*%B + z%*%matrix(1,ncol = J, nrow = 1)
    lls[iter + 1] <- log_likelihood_wide(Y = Y,
                                         rect_weights = rect_weights,
                                         X = X,
                                         B = B,
                                         z = z)


    if((abs(lls[iter + 1] - lls[iter]) < tolerance)&
       (max(lls[1:iter] - lls[iter + 1]) < tolerance)
    ){
      converged <- TRUE
    }

    if(iter>=maxit){
      converged <- TRUE
    }



    iter <- iter + 1
    # par(mfrow =  c(2,1))
    # plot(asinh(lls[1:(iter)]),type = "l")
    # plot(asinh(lls[1:(iter)] - lls[iter]),type = "l")
  }

  if(verbose){
    message(paste("Iteration ", iter,"; log likelihood ", lls[iter],
                  sep ="", collapse = ""))
  }

  if(return_a_lot){
    return(list("B" = B,
                "z" = z,
                "ll" = lls[iter - 1],
                "Y" = Y,
                "X" = X,
                "tolerance" = tolerance,
                "maxit" = maxit,
                "constraint_fn" = constraint_fn,
                "maxit_glm" = maxit_glm,
                "method" = method,
                "weights" = rect_weights,
                "reweight" = reweight))
  } else{
    return(
      list("B" = B,
           "z" = z,
           "ll" = lls[iter - 1],
           "method" = method,
           "weights" = rect_weights,
           "reweight" = reweight)
    )

  }
}
