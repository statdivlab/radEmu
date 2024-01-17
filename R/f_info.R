#' @importFrom methods as

f_info <- function(Y,
                   B_cup,
                   B,
                   X,
                   X_cup,
                   omit = NULL,
                   compute_together = FALSE,
                   return_separate = FALSE){
  n <- nrow(Y)
  J <- ncol(Y)
  p <- ncol(X)

  if(is.null(omit)){
  f_info <- matrix(0,ncol = p*J,nrow = p*J)
  } else{
    if(compute_together){
    f_info <- matrix(0,ncol = p*(J - 1), nrow = p*(J - 1))
    }
  }
  z <- update_z(Matrix::Matrix(Y),Matrix::Matrix(X),B)



  if(compute_together){
    ident_mat <- Matrix::Diagonal(x = rep(1,J))
    one_mat <- Matrix::Matrix(1,nrow = J,ncol=1)
  for(i in 1:n){
    print(i)
    which_rows <- 1:J + (i - 1)*J
    mu_i <- exp(X_cup[which_rows,]%*% B_cup + z[i])
    p_i <- mu_i/sum(mu_i)

    M_iX_i <- (ident_mat - Matrix::tcrossprod(one_mat,p_i))%*%X_cup[which_rows,]


    f_info <- f_info + Matrix::crossprod(M_iX_i,Matrix::Diagonal(x= as.numeric(mu_i)))%*%M_iX_i
  }
  } else{
    mus <- exp(X_cup%*%B_cup + rep(z,each = J))
    mu_diag <- Matrix::Diagonal(x = mus@x)

    ps <- lapply(1:n,
                 function(i){
                   p_i <- mus@x[(i-1)*J +1:J];
                 return(p_i/sum(p_i))})
    ps_mat <- Matrix::sparseMatrix(i = 1:(n*J),
                                   j = rep(1:n,each = J),
                                   x = do.call(c,ps))
    X_cup_p <- Matrix::crossprod(X_cup,ps_mat)
    X_cup_p <- as.matrix(X_cup_p)
    X_cup_p <- methods::as(X_cup_p,"sparseMatrix")
    sum_diag <- Matrix::Diagonal(x = rowSums(Y))

    f_info <- Matrix::crossprod(X_cup,mu_diag)%*%X_cup -
      X_cup_p%*%
      Matrix::tcrossprod(sum_diag,X_cup_p)



    }



  return(f_info)


}
