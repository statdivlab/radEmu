null_repar_ll <- function(x,js,B,z,p,Y,X,j_constr,k_constr,constraint_fn){
  # print(get_ll(Y = Y, X = X, B = B, z = z))
  # 
  # njs <- length(js)
  # 
  # for(jind in 1:njs){
  #   temp_j <- js[jind]
  #   B[,temp_j] <- B[,temp_j] + x[1:p + (jind - 1)*p]
  # }
  # 
  # B[-k_constr,j_constr] <- B[-k_constr,j_constr] + x[(njs*p) +( 1:(p - 1))]
  # B[k_constr,j_constr] <- constraint_fn(B[k_constr,-j_constr])
  # return(get_ll(Y = Y,X = X,B = B,z = z))
  Bjs <- B[,c(js,j_constr),drop = FALSE]
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
  # 
  return(sum(Y[,c(js,j_constr)]*log_means_js - exp(log_means_js)))
  
}