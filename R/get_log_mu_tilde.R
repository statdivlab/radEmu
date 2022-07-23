

get_log_mu_tilde <- function(X_tilde,B_tilde,n,J,A = NULL){

log_mu <- X_tilde%*%B_tilde

if(is.null(A)){
mini_A <- matrix(1,nrow = J, ncol = J)
A <- Matrix::bdiag(lapply(1:n, function(i) mini_A))
}

log_mu_tilde <- log_mu - log(A%*%exp(log_mu))

return(log_mu_tilde)
}
