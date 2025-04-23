
get_constrained_gr <- function(Y,X,B,z,js,j_constr,k_constr,j_ref,
                               constraint_grad_fn,gr_only = FALSE){
  Bjs <- B[,c(js,j_constr)]
  Yjs <- Y[,c(js,j_constr)]
  log_means_js <- X%*%Bjs
  for(i in 1:nrow(Y)){
    log_means_js[i,] <- log_means_js[i,] + z[i]
  }
  njs <- length(js)
  
  gr <- lapply(1:(njs + 1),
               function(jind){
                 Matrix::crossprod(X,
                                   Yjs[,jind,drop = FALSE] - 
                                     exp(log_means_js[,jind,drop = FALSE]))
               }
               )
  cg <- constraint_grad_vec(constraint_grad_fn,
                            js_used = js,
                            Bk_constr = B[k_constr,],
                            j_constr = j_constr,
                            p = ncol(X))
  
  cg_gr_multiplier <- gr[[length(gr)]][k_constr]
  
  for(jind in 1:njs){
    gr[[jind]][k_constr] <-  gr[[jind]][k_constr] +
      cg[jind]*cg_gr_multiplier

  }
  gr[[length(gr)]] <- gr[[length(gr)]][-k_constr]
  gr <- do.call(c,gr)
  
  if(gr_only){
    return(gr)
  }
  # 
  # r_info <- restr_info(log_means_js,js,j_constr,k_constr,j_ref,J)
  # 
  # info_j_constr_k_constr <- sum((X[k_constr,])^2*exp(log_means_js[,njs + 1]))
  # 
  # # for(jind in 1:njs){
  # # 
  # #   to_add <- -numDeriv::grad(function(x){
  # #     temp_Bkconstr <- B[k_constr,];
  # #     temp_Bkconstr[js[jind]] <- x;
  # #     temp_Bkconstr <- temp_Bkconstr[-j_constr];
  # #     ind <- js[ind] - as.numeric(js[ind] > j_constr);
  # #     return(constraint_grad_fn(temp_Bkconstr)[ind])
  # #   },      B[k_constr,js[jind]])
  # #   to_add <- to_add* sum(X[k_constr,]*(Yjs[,njs + 1] - exp(log_means_js[,njs + 1])))
  # #   r_info[(jind - 1)*p + k_constr,(jind - 1)*p + k_constr] <-
  # #     r_info[(jind - 1)*p + k_constr,(jind - 1)*p + k_constr] + to_add
  # #   # print(to_add)
  # # }
  # # 
  # # r_info_inv <- Matrix::solve(r_info)
  # # 
  # # 
  # # 
  # # I_inv_nab_l <- r_info_inv%*%matrix(gr,ncol = 1)
  # # 
  # u <- cg*info_j_constr_k_constr
  # u <- lapply(1:njs,function(jind){
  #   x <- rep(0,p);
  #   x[k_constr] <- u[jind];
  #   return(x)
  # })
  # u <- matrix(c(do.call(c,u),0),ncol = 1)
  # 
  # cg <- constraint_grad_fn(B[k_constr,c(js,j_ref)])[1:(J - 2)]
  # UUT <- 0*r_info
  # for(jind in 1:njs){
  #   to_add <- sapply(1:njs,function(ind) -numDeriv::grad(function(x){
  #         temp_Bkconstr <- B[k_constr,];
  #         temp_Bkconstr[js[jind]] <- x;
  #         temp_Bkconstr <- temp_Bkconstr[-j_constr];
  #         ind <- js[ind] - as.numeric(js[ind] > j_constr);
  #         return(constraint_grad_fn(temp_Bkconstr)[ind])
  #       },      B[k_constr,js[jind]]))*sum(X[,k_constr]*(Yjs[,njs + 1]- exp(log_means_js[,njs + 1])))
  #   for(jind2 in 1:njs){
  #     UUT[p*(jind2 - 1) + k_constr,p*(jind - 1) + k_constr] <- 
  #       UUT[p*(jind - 1) + k_constr,p*(jind2 - 1) + k_constr] <- 
  #       to_add[jind2]
  #   }
  # }
  # # 
  # a <- r_info_inv%*%u
  # to_subtract <-
  #   a*(sum(a*gr))/(1 + sum(u*a))
  # 
  # return(list("update" = I_inv_nab_l - to_subtract,"grad" = gr,"info" = r_info + Matrix::tcrossprod(matrix(u,ncol = 1))))
  # r_info <- r_info + Matrix::tcrossprod(u) + UUT
  # return(list("update" = qr.solve(r_info,matrix(gr,ncol = 1)),
  #             "grad" = gr))
}