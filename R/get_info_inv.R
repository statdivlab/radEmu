

get_info_inv <- function(X_tilde_repar,
                         S,
                         J,
                         n,
                         glm_weights){
  V <- Matrix::Diagonal(x = glm_weights)

  Amat <-  Matrix::crossprod(V%*%X_tilde_repar,X_tilde_repar)

A_inv <- Matrix::solve(
 Amat
)

  Bmat <- Matrix::crossprod(X_tilde_repar,V%*%S)
  A_inv_B <- Matrix::solve(Amat,Bmat)

  C_A_inv_B <- Matrix::crossprod(Bmat,A_inv_B)
  D <- Matrix::crossprod(V%*%S,S)
  D_minus_C_A_inv_B_inv <- Matrix::solve(D- C_A_inv_B)

  part11 <- A_inv -
                   A_inv%*%Bmat%*%(Matrix::tcrossprod(D_minus_C_A_inv_B_inv,
                                                  Bmat))%*%A_inv
  part21 <- -Matrix::tcrossprod(D_minus_C_A_inv_B_inv,Bmat)%*%A_inv
  part12 <- Matrix::t(part21)
  part22 <- D_minus_C_A_inv_B_inv

  return(cbind(rbind(part11,part21),
               rbind(part12,part22)))

}
