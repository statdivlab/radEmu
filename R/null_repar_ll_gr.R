null_repar_ll_gr <- function(x,js,B,z,p,Y,X,j_constr,k_constr,constraint_fn,
                      constraint_grad_fn,
                      return_hess = FALSE){
  Bjs <- B[,c(js,j_constr)]
  njs <- length(js)
  for(jind in 1:njs){
    Bjs[,jind] <- Bjs[,jind] + x[1:p + (jind -1)*p]
  }
  Bjs[-k_constr,njs + 1] <- Bjs[-k_constr,njs + 1] + x[p*njs + 1:(p-1)]
  Bjs[k_constr,njs + 1] <- constraint_fn(c(B[k_constr,-c(js,j_constr)],
                                           Bjs[k_constr,1:njs]))
  # 
  log_means_js <- X%*%Bjs
  # 
  for(i in 1:nrow(Y)){
    log_means_js[i,] <- log_means_js[i,] + z[i]
  }
  
  Yjs <- Y[,c(js,j_constr)]
  # 
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
  
  if(return_hess){
    info_diags <- lapply(1:njs,
                         function(jind) Matrix::crossprod(X,
                                                          diag(exp(log_means_js[,jind])))%*%X)
    
    
    info_diags <- c(info_diags,
                    Matrix::crossprod(X[,-k_constr,drop = FALSE],
                                      diag(exp(log_means_js[,njs + 1])))%*%X[,-k_constr,
                                                                             drop = FALSE])
    info <- Matrix::bdiag(info_diags)
    info_inv <- Matrix::bdiag(lapply(info_diags,qr.solve))
    
    return(list(gr = gr,
                info = info,
                info_inv = info_inv,
                cg = cg,
                cg_info_mult = Matrix::crossprod(X[,k_constr,drop = FALSE],
                                                 diag(exp(log_means_js[,njs + 1])))%*%
                  X[,k_constr,drop = FALSE]))
    
  } else{
    
    
    return(gr)}
}