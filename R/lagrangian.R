

lagrangian <- function(X,Y,B,lambda_ks,lambda_kj, null_k, null_j){
  dBcup = dpll_dB_cup(X,Y,B)

  p <- ncol(X)
  if(length(lambda_ks)!= p){
    stop("lambda_ks must have number of entries equal to ncol(X)")
  }

  null_index <- p*(null_j - 1) + null_k

  lagrangian_vec = c(as.numeric(dBcup), rep(0,p + 1))

  lagrangian_vec <- lagrangian_vec + c(rep(lambda_ks,J),rowSums(B),B[null_k,null_j])

  lagrangian_vec[null_index] <-  lagrangian_vec[null_index] + lambda_kj

  return(lagrangian_vec)
}
