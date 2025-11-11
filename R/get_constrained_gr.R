get_constrained_gr <- function(Y, X, B, z, js, j_constr, k_constr, j_ref,
                               constraint_fn, constraint_grad_fn, 
                               gr_only = FALSE, ref_set = NULL){
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
  cg <- constraint_grad_vec(constraint_fn,
                            constraint_grad_fn,
                            js_used = js,
                            Bk_constr = B[k_constr,],
                            j_constr = j_constr,
                            p = ncol(X),
                            ref_set = ref_set)
  
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
}