

get_mu_tilde <- function(rho,X,Y){
  J <- ncol(rho)
  denoms <- sapply(1:n,function(i) logsum::sum_of_logs(
    sapply(1:J, function(j) X[i,,drop=FALSE]%*%rho[,j,drop=FALSE])))
return(exp(X%*%rho - matrix(denoms,ncol = 1)%*%matrix(1,nrow = 1, ncol = J)))
}
