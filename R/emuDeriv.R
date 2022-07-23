

emuDeriv <- function(D_tilde,
                     theta,
                     Y_tilde,
                     weights){


return(Matrix::crossprod(Y_tilde - exp(D_tilde%*%theta),
                  Matrix::Diagonal(x = weights)%*%D_tilde))
}
