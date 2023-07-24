lagrangian_full_model <- function(X,Y,B,lambda_ks){
  dBcup = dpll_dB_cup(X,Y,B)

  p <- ncol(X)
  J <- ncol(Y)
  if(length(lambda_ks)!= p){
    stop("lambda_ks must have number of entries equal to ncol(X)")
  }


  lagrangian_vec = c(as.numeric(dBcup), rep(0,p ))

  lagrangian_vec <- lagrangian_vec + c(rep(lambda_ks,J),rowSums(B))

  return(lagrangian_vec)
}
