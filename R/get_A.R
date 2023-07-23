

get_A <- function(B_cup,
                  X_cup,
                  Y){
  n <- nrow(Y)

  X_cup_i <- lapply(1:n,
                    function(i)
                      X_cup[(i - 1)*J + 1:J,])

  z_i <- lapply(1:n,
                function(i)
                 rep(log(sum(Y[i,])) - log(sum(exp( X_cup_i[[i]]%*%B_cup))),
                     J))

mu_i <- lapply(1:n,
               function(i)
                 as.numeric(exp(X_cup_i[[i]]%*%B_cup + z_i[[i]])))

Dlog_mu_i_dB_cup <-
  lapply(1:n,function(i)
    X_cup_i[[i]]
    -matrix(1,ncol = 1,nrow = J)%*%matrix(Matrix::colSums(

      Matrix::Diagonal(x = mu_i[[i]]/sum(mu_i[[i]]))%*%
                                                    X_cup_i[[i]]),nrow = 1))



 A_parts <-
   lapply(1:n, function(i)
     Matrix::crossprod(Dlog_mu_i_dB_cup[[i]],
                       Matrix::Diagonal(x= mu_i[[i]])) %*%
                               Dlog_mu_i_dB_cup[[i]]
   )

A <- A_parts[[1]]
for(i in 2:n){
  A <- A + A_parts[[i]]
}

return(list("A" = A,
            "Dlog_mu_i_dB_cup" = Dlog_mu_i_dB_cup))

}
