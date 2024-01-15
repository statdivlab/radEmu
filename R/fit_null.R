#' Fit model under null
#'
#'
fit_null <- function(B,
                     #B (MPLE)
                     Y, #Y (with augmentations)
                     X, #design matrix
                     X_cup = NULL,
                     k_constr, #row index of B to constrain
                     j_constr, #col index of B to constrain
                     j_ref,
                     constraint_fn, #constraint function
                     constraint_grad_fn, #gradient of constraint fn\
                     rho_init = 1,
                     tau = 1.2,
                     kappa = 0.8,
                     B_tol = 1e-2,
                     inner_tol = 0.01,
                     constraint_tol = 1e-4,
                     max_step = 5,
                     c1 = 1e-4,
                     maxit = 1000,
                     inner_maxit = 25,
                     n_final_inner_it = 25,
                     verbose = FALSE,
                     trackB = FALSE
){

  J <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)

  #store cols of B to update in a vector
  j_to_update <- (1:J)[(1:J) != j_ref]
  #in vector format, which indexes to remove (bc B^j_constr = 0 id. constr.)
  indexes_to_remove <- (j_ref - 1)*p + 1:p


  #change id. constr. by subtracting approp. col of B from others
  for(k in 1:p){
    B[k,] <- B[k,] - B[k,j_ref]
  }

  #how large is norm of score at initiation?
  z <- update_z(X = X, Y = Y, B = B)

  log_mean <- X%*%B +
    matrix(z,ncol = 1)%*%matrix(1,ncol = J, nrow = 1)

  # deriv <- do.call(c,lapply(1:J,
  #                           function(j) crossprod(X,Y[,j,drop = FALSE] - exp(log_mean[,j,drop = FALSE]))))
  #
  # inner_tol <-  inner_tol*sqrt((1/(n*J))*sum(deriv^2))

  # message("Setting tolerance for inner iterations to ", signif(inner_tol,2),".")

  if(trackB){
    Bs <- data.frame(k = rep(1:p,each = J),
                     j = rep(1:J,p),
                     B = NA)
  } else{
    Bs <- NULL
  }

  #get X_cup for later use
  if(is.null(X_cup)){
    X_cup = X_cup_from_X(X,J)
  }

  #set  <- ation to zero
  iter <- 0

  #compute gap
  gap <- constraint_fn(B[k_constr,]) - B[k_constr,j_constr]
  init_gap <- abs(gap)

  #set rho equal to initial value
  rho <- rho_init

  #initiate u
  u <- rho*gap

  #deriv of (poisson) ll with z fixed

  # dldB <- lapply(1:J,function(j) if(j != j_constr){
  #   means <- exp(X%*%B[,j] + z)
  #   return(crossprod(X,Y[,j] - means))} else{NULL}
  # )
  #
  # dldB <- do.call(cbind,dldB)
  #
  # cg <- constraint_grad_fn(B[k_constr,])
  # means <- exp(X%*%B[,j,drop = FALSE] + z)
  # dldB[k_constr,] <-   dldB[k_constr,] +
  #   (u + rho*gap)*cg[-j_constr]
  #
  # # Q <- sqrt((1/(n*J))*sum(dldB^2))
  # Q <- Inf

  #
  # entries_for_g <- rep((1:p) == k_constr,J)[-indexes_to_remove]
  # dgdB <- 0*dldB
  # dgdB[entries_for_g,] <- constraint_grad_fn(B[k_constr,])[-j_constr]
  #
  # lambda_hat <- -as.numeric(crossprod(dldB,dgdB)/crossprod(dgdB,dgdB))
  #
  # dlag_lamhat_dB <- dldB + lambda_hat*dgdB
  #
  # Q <- sum(dlag_lamhat_dB^2)

  if(trackB){
    Bs$iter <- 0
    Bs$inner_iter <- 0
    Bs$rho <- rho
    Bs$u <- u
    Bs$Q <- NA
    Bs$gap <- gap
    Bs$score_stat <- 0
    Bs$score_err <- 0
    Bs$dlag_dB <- NA
    lilBs <- Bs

    for(k in 1:p){
      Bs$B[1:J + J*(k - 1)] <- B[k,]
      # print(unique(Bs$k[1:J + J*(k - 1)]))
    }}

  loop_j_to_update <- j_to_update

  # proj_score <- Q <- Inf
  use_max_inner_it <- FALSE
  B_diff <- Inf

  while((abs(gap) > constraint_tol | B_diff> B_tol #proj_score > B_tol #| Q > inner_tol
         ) & iter <= maxit){
    #increment iteration
    iter <- iter + 1
    #initiate internal iteration
    inner_iter <- 0

    #evaluate augmented Lagrangian
    log_means <- do.call(cbind,lapply(1:J,function(j) X%*%B[,j] + z))

#     message("Have to update how lagrangian and its derivs are computed to reflect
# better ref taxon")
    curr_lag_val <- sum(-log_means*Y + exp(log_means)) + u*gap + (rho/2)*(gap^2)

    old_gap <- gap
    old_B <- B
    inner_diff <- Inf
    while(((inner_diff > inner_tol|use_max_inner_it) & inner_iter <= inner_maxit)| inner_iter == 0){
      inner_old_B <- B
      # column-wise optimization:
        # dldB <- lapply(1:J,function(j) if(j != j_constr){
        #   means <- exp(X%*%B[,j] + z)
        #   return(-crossprod(X,Y[,j] - means))} else{NULL}
        # )
        #
        # dldB <- do.call(cbind,dldB)
        #
        # cg <- constraint_grad_fn(B[k_constr,])
        # means <- exp(X%*%B[,j,drop = FALSE] + z)
        # dldB[k_constr,] <-   dldB[k_constr,] +
        #   (u + rho*gap)*cg[-j_constr]
        # cs <- colSums(dldB^2)
        # loop_j_to_update <- which(cs>=quantile(cs,0.95))
        # loop_j_to_update[loop_j_to_update >= j_constr] <- loop_j_to_update[loop_j_to_update >= j_constr] + 1
        # for(j in loop_j_to_update){
        #
        #   # print(j)
        #   #compute derivatives
        #   # cg <- constraint_grad_fn(B[k_constr,])
        #   # ch <- constraint_hess_fn(B[k_constr,],j,j)
        #   # means <- exp(X%*%B[,j,drop = FALSE] + z)
        #   # lag_grad_in_Bj <- -crossprod(X,Y[,j] - means)
        #   # lag_grad_in_Bj[k_constr,] <- lag_grad_in_Bj[k_constr,] +
        #   #   (u + rho*gap)*cg[j]
        #   #
        #   # # num_grad_Bj <- numDeriv::grad(function(x){
        #   # #   log_means <- X%*%matrix(x,ncol = 1) + z;
        #   # #   Bkconstr <- B[k_constr,];
        #   # #   Bkconstr[j] <- x[k_constr];
        #   # #   gp <- constraint_fn(Bkconstr);
        #   # #   return(sum(-Y[,j]*log_means + exp(log_means)) + u*gp + (rho/2)*gp^2)
        #   # # },
        #   # # B[,j])
        #   #
        #   # lag_hess_in_Bj <- crossprod(X,diag(as.numeric(means)))%*%X
        #   # lag_hess_in_Bj[k_constr,k_constr] <-
        #   #   lag_hess_in_Bj[k_constr,k_constr] + rho*(cg[j])^2 +
        #   #   ch*(u + rho*gap)
        #   #
        #   #
        #   # update_direction <- try(-qr.solve(lag_hess_in_Bj,
        #   #                                   lag_grad_in_Bj,
        #   #                                   tol = 1e-20),silent = TRUE)
        #   # if(inherits(update_direction,"try-error")){
        #   #   update_direction <- -lag_grad_in_Bj
        #   #   if(verbose){
        #   #     message("Using gradient update")}
        #   # }
        #   #
        #   # armijo_term <- c1*sum(lag_grad_in_Bj*update_direction)
        #   #
        #   # if(armijo_term >0){
        #   #   update_direction <- -lag_grad_in_Bj
        #   #   if(verbose){
        #   #     message("Using gradient update")}
        #   #   armijo_term <-  c1*sum(lag_grad_in_Bj*update_direction)
        #   # }
        #   #
        #   # # if(max(abs(update_direction))>.1){
        #   # #   update_direction <- .1*update_direction/max(abs(update_direction))
        #   # #   if(verbose){
        #   # #   message("Update step has l_inf norm >.1; rescaling so norm is .1")}
        #   # #   armijo_term <-  c1*sum(lag_grad_in_Bj*update_direction)
        #   # # }
        #   # accept_step <- FALSE
        #   #
        #   # stepsize <- 0.5
        #   # while((!accept_step)){
        #   #   prop_B <- B
        #   #   prop_B[,j] <- prop_B[,j] + stepsize*update_direction
        #   #   prop_gap <- constraint_fn(prop_B[k_constr,])
        #   #   prop_log_means <- do.call(cbind,lapply(1:J,function(j) X%*%prop_B[,j] + z))
        #   #   prop_lag_val <- sum(-prop_log_means*Y + exp(prop_log_means)) +
        #   #     u*prop_gap + (rho/2)*(prop_gap^2)
        #   #   accept_step <- prop_lag_val <= curr_lag_val + stepsize*armijo_term
        #   #   # accept_step <- TRUE
        #   #   stepsize <- stepsize/2
        #   #
        #   #   if(stepsize<1e-15){
        #   #     prop_B <- B
        #   #     accept_step <- TRUE
        #   #   }
        #   # }
        #   # B[,j] <- prop_B[,j]
        #   # gap <- prop_gap
        #   # print(j)
        #
        #   update <- micro_fisher_null(X = X,
        #                               Y = Y,
        #                               B = B,
        #                               j = j,
        #                               u = u,
        #                               rho = rho,
        #                               k_constr = k_constr,
        #                               constraint_fn = constraint_fn,
        #                               constraint_grad_fn = constraint_grad_fn,
        #                               constraint_hess_fn = constraint_hess_fn,
        #                               z = z,
        #                               stepsize = 1,
        #                               c1 = c1)
        #
        #   B[,j] <- B[,j] + update$update
        #   gap <- update$updated_gap
        #
        #   # message("Step for column ", j," accepted with stepsize ",stepsize*2)
        # }
        #
        # z <- update_z(Y,X,B)
#
#       message("Have to update how lagrangian and its derivs are computed to reflect
# better ref taxon")
      update <- macro_fisher_null(X = X,
                                  Y = Y,
                                  B = B,
                                  z = z,
                                  J = J,
                                  p = p,
                                  k_constr = k_constr,
                                  j_constr = j_constr,
                                  j_ref = j_ref,
                                  rho = rho,
                                  u = u,
                                  max_step = max_step,
                                  constraint_fn = constraint_fn,
                                  constraint_grad_fn = constraint_grad_fn,
                                  verbose = verbose, 
                                  c1 = 1e-4)
      B <- B + update$update
      gap <- update$gap
      z <- update_z(Y,X,B)


      if(trackB){
        temp_Bs <- lilBs
        for(k in 1:p){

          temp_Bs$B[1:J + J*(k - 1)] <- B[k,]

        }
        temp_Bs$iter <- iter
        temp_Bs$inner_iter <- inner_iter
        temp_Bs$rho <- rho
        temp_Bs$u <- u
        temp_Bs$Q <- Q
        temp_Bs$gap <- gap
        if(iter>1|inner_iter>1){
        jcounter <- 1
        for(j in 1:J){
          if(j == j_constr){
            temp_Bs$dlag_dB[temp_Bs$j ==j] <- 0
          } else{
            temp_Bs$dlag_dB[temp_Bs$j ==j] <- dlag_lamhat_dB[(jcounter - 1)*p + 1:p,]
            jcounter <- jcounter + 1
          }
        }}


        Bs <- rbind(Bs,temp_Bs)



        # print(Q)

      }

      inner_diff <- max(abs(B - inner_old_B))

      #increment inner loop iteration counter
      inner_iter <- inner_iter + 1
      # print(inner_iter)


#
#       message("Have to update how lagrangian and its derivs are computed to reflect
# better ref taxon")
      #evaluate gradient of augmented Lagrangian and also projected score
      # dldB <- lapply(1:J,function(j) if(j != j_ref){
      #   means <- exp(X%*%B[,j] + z)
      #   return(crossprod(X,Y[,j] - means))} else{NULL}
      # )
      #
      # dldB_for_hat <- do.call(rbind,dldB)
      # dldB <- do.call(cbind,dldB)
      #
      # entries_for_g <- rep((1:p) == k_constr,J)[-indexes_to_remove]
      # dgdB <- 0*dldB_for_hat
      # cg <- constraint_grad_fn(B[k_constr,])
      # cg[j_constr] <- cg[j_constr] - 1
      # dgdB[entries_for_g,] <- cg[-j_ref]
      #
      # lambda_hat <- -as.numeric(crossprod(dldB_for_hat,dgdB)/crossprod(dgdB,dgdB))
      #
      #
      # dlag_lamhat_dB <- dldB_for_hat + lambda_hat*dgdB
      #
      # dldB[k_constr,] <-   dldB[k_constr,] +
      #   (u + rho*gap)*cg[-j_ref]
      #
      # Q <- min(sqrt((1/(n*J))*sum(dlag_lamhat_dB^2)),
      #          sqrt((1/(n*J))*sum(dldB^2)))


    }

    B_diff <- max(abs(B - old_B))


    if(verbose){
      message("Max absolute difference in B since last augmented Lagrangian outer step: ", signif(B_diff,3))
    }


    if( abs(gap)>constraint_tol){
      u <- u + rho*gap
      if(abs(gap/old_gap) > kappa){
        rho <- rho*tau
      } else{
        # rho <- rho*(tau^(-0.25))
      }
    }




      if(verbose){
        message("Estimate deviates from feasibility by ",signif(gap,3))
        message("Parameter u set to ", signif(u,3))
        message("Parameter rho set to ", signif(rho,3))
    }



#     message("Have to update how lagrangian and its derivs are computed to reflect
# better ref taxon")
    # dldB <- lapply(1:J,function(j) if(j != j_constr){
    #   means <- exp(X%*%B[,j] + z)
    #   return(crossprod(X,Y[,j] - means))} else{NULL}
    # )
    #
    # dldB_for_hat <- do.call(rbind,dldB)
    # dldB <- do.call(cbind,dldB)
    #
    # entries_for_g <- rep((1:p) == k_constr,J)[-indexes_to_remove]
    # dgdB <- 0*dldB_for_hat
    # dgdB[entries_for_g,] <- constraint_grad_fn(B[k_constr,])[-j_constr]
    #
    # lambda_hat <- -as.numeric(crossprod(dldB_for_hat,dgdB)/crossprod(dgdB,dgdB))
    #
    #
    # dlag_lamhat_dB <- dldB_for_hat + lambda_hat*dgdB
    #
    # Q <- min(sqrt((1/(n*J))*sum(dlag_lamhat_dB^2)),
    #          sqrt((1/(n*J))*sum(dldB^2)))
    if(verbose){

      # message("Convergence criterion equal to ", signif(Q,3))
    }


    if(abs(gap) < constraint_tol){
#       message("Have to update how lagrangian and its derivs are computed to reflect
# better ref taxon")
    # check <- get_score_stat(Y = Y,
    #                X_cup = X_cup,
    #                X = X,
    #                B = B,
    #                k_constr = k_constr,
    #                j_constr = j_constr,
    #                j_ref = j_ref,
    #                constraint_grad_fn = constraint_grad_fn,
    #                indexes_to_remove = indexes_to_remove,
    #                J = J,
    #                n = n,
    #                p = p,
    #                check_influence = TRUE,
    #                I = I,
    #                Dy = Dy)

    use_max_inner_it <- TRUE



    # message("Score stat: ",signif(check$score_stat*(n/(n - 1)),3))
    # message("Influence stat: ",signif(check$influence_stat,3))
    # proj_score <- as.numeric(as.matrix(check$influence_stat))
    } else{
      use_max_inner_it <- FALSE
    }


}

  # print(Q)
  # print(inner_tol)

  return(list("B" = B,
              "k_constr" = k_constr,
              "j_constr" = j_constr,
              # "score_stats" = list(score_stats),
              # "proj_score" = proj_score,
              "niter" = iter,
              "gap" = gap,
              "u" = u,
              "rho" = rho,
              "Bs" = Bs))




}

# column-wise optimization:
# dldB <- lapply(1:J,function(j) if(j != j_constr){
#   means <- exp(X%*%B[,j] + z)
#   return(-crossprod(X,Y[,j] - means))} else{NULL}
# )
#
# dldB <- do.call(cbind,dldB)
#
# cg <- constraint_grad_fn(B[k_constr,])
# means <- exp(X%*%B[,j,drop = FALSE] + z)
# dldB[k_constr,] <-   dldB[k_constr,] +
#   (u + rho*gap)*cg[-j_constr]
# cs <- colSums(dldB^2)
# loop_j_to_update <- which(cs>=quantile(cs,0.95))
# loop_j_to_update[loop_j_to_update >= j_constr] <- loop_j_to_update[loop_j_to_update >= j_constr] + 1
# for(j in loop_j_to_update){
#
#   # print(j)
#   #compute derivatives
#   # cg <- constraint_grad_fn(B[k_constr,])
#   # ch <- constraint_hess_fn(B[k_constr,],j,j)
#   # means <- exp(X%*%B[,j,drop = FALSE] + z)
#   # lag_grad_in_Bj <- -crossprod(X,Y[,j] - means)
#   # lag_grad_in_Bj[k_constr,] <- lag_grad_in_Bj[k_constr,] +
#   #   (u + rho*gap)*cg[j]
#   #
#   # # num_grad_Bj <- numDeriv::grad(function(x){
#   # #   log_means <- X%*%matrix(x,ncol = 1) + z;
#   # #   Bkconstr <- B[k_constr,];
#   # #   Bkconstr[j] <- x[k_constr];
#   # #   gp <- constraint_fn(Bkconstr);
#   # #   return(sum(-Y[,j]*log_means + exp(log_means)) + u*gp + (rho/2)*gp^2)
#   # # },
#   # # B[,j])
#   #
#   # lag_hess_in_Bj <- crossprod(X,diag(as.numeric(means)))%*%X
#   # lag_hess_in_Bj[k_constr,k_constr] <-
#   #   lag_hess_in_Bj[k_constr,k_constr] + rho*(cg[j])^2 +
#   #   ch*(u + rho*gap)
#   #
#   #
#   # update_direction <- try(-qr.solve(lag_hess_in_Bj,
#   #                                   lag_grad_in_Bj,
#   #                                   tol = 1e-20),silent = TRUE)
#   # if(inherits(update_direction,"try-error")){
#   #   update_direction <- -lag_grad_in_Bj
#   #   if(verbose){
#   #     message("Using gradient update")}
#   # }
#   #
#   # armijo_term <- c1*sum(lag_grad_in_Bj*update_direction)
#   #
#   # if(armijo_term >0){
#   #   update_direction <- -lag_grad_in_Bj
#   #   if(verbose){
#   #     message("Using gradient update")}
#   #   armijo_term <-  c1*sum(lag_grad_in_Bj*update_direction)
#   # }
#   #
#   # # if(max(abs(update_direction))>.1){
#   # #   update_direction <- .1*update_direction/max(abs(update_direction))
#   # #   if(verbose){
#   # #   message("Update step has l_inf norm >.1; rescaling so norm is .1")}
#   # #   armijo_term <-  c1*sum(lag_grad_in_Bj*update_direction)
#   # # }
#   # accept_step <- FALSE
#   #
#   # stepsize <- 0.5
#   # while((!accept_step)){
#   #   prop_B <- B
#   #   prop_B[,j] <- prop_B[,j] + stepsize*update_direction
#   #   prop_gap <- constraint_fn(prop_B[k_constr,])
#   #   prop_log_means <- do.call(cbind,lapply(1:J,function(j) X%*%prop_B[,j] + z))
#   #   prop_lag_val <- sum(-prop_log_means*Y + exp(prop_log_means)) +
#   #     u*prop_gap + (rho/2)*(prop_gap^2)
#   #   accept_step <- prop_lag_val <= curr_lag_val + stepsize*armijo_term
#   #   # accept_step <- TRUE
#   #   stepsize <- stepsize/2
#   #
#   #   if(stepsize<1e-15){
#   #     prop_B <- B
#   #     accept_step <- TRUE
#   #   }
#   # }
#   # B[,j] <- prop_B[,j]
#   # gap <- prop_gap
#   # print(j)
#
#   update <- micro_fisher_null(X = X,
#                               Y = Y,
#                               B = B,
#                               j = j,
#                               u = u,
#                               rho = rho,
#                               k_constr = k_constr,
#                               constraint_fn = constraint_fn,
#                               constraint_grad_fn = constraint_grad_fn,
#                               constraint_hess_fn = constraint_hess_fn,
#                               z = z,
#                               stepsize = 1,
#                               c1 = c1)
#
#   B[,j] <- B[,j] + update$update
#   gap <- update$updated_gap
#
#   # message("Step for column ", j," accepted with stepsize ",stepsize*2)
# }
#
# z <- update_z(Y,X,B)
#
