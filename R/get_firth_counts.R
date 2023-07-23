get_firth_counts <- function(Y,
                             B_cup,
                             B,
                             X,
                             X_cup){
  n <- nrow(Y)
  p <- ncol(X)
  J <- ncol(Y)

 W <- f_info_ingredients(Y,B_cup,B,X,X_cup)

 W_eigen <- eigen(W)
 W_sqrt <- W_eigen$vectors%*%diag(sqrt(pmax(W_eigen$values,0)))%*%t(W_eigen$vectors)

 n_to_keep <- p*(J - 1)
 n_remaining <- p
 info <- Matrix::t(X_cup)%*%W%*%X_cup
 info_eigen <- eigen(info)
 info_neg_sqrt <- info_eigen$vectors %*%
   diag(c(1/sqrt(info_eigen$values[1:n_to_keep]),rep(0,n_remaining))) %*%
   t(info_eigen$vectors)

 scaled_X_cup <- W_sqrt%*%X_cup%*%info_neg_sqrt

 return(0.5*Matrix::rowSums(scaled_X_cup^2))

}
