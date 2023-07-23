

constr_lagrangian <- function(Y,
                              X,
                              B,
                              constraint_fn,
                              constraint_grad_fn,
                              k_constr,
                              j_constr
){

  ll_grad <- dpll_dB_cup(X,Y,B)

  e_k_constr_j_constr <- 0*B
  e_k_constr_j_constr[k_constr,j_constr] <- 1
  e_k_constr_j_constr <- B_cup_from_B(  e_k_constr_j_constr)

  constraint_grad <- 0*B
  constraint_grad[k_constr,-j_constr] <- constraint_grad_fn(B[k_constr,-j_constr])
  constraint_grad <- B_cup_from_B(constraint_grad)

  lambda_numerator <- -1*as.numeric(
    as.matrix(
      Matrix::crossprod( e_k_constr_j_constr - constraint_grad,ll_grad)
      )
    )
  lambda_denominator <- as.numeric(as.matrix(Matrix::crossprod(e_k_constr_j_constr  - constraint_grad)))
  lambda <- lambda_numerator/lambda_denominator


  lagrangian_deriv <- ll_grad + lambda*(e_k_constr_j_constr - constraint_grad)
  return(lagrangian_deriv)
}


# lambda_func <- function(lamb){
#   return(sum((ll_grad + lamb*(e_k_constr_j_constr - constraint_grad))^2))
# }
#
# lambdas <- seq(-1000,2000,100)
#
# norms <- sapply(lambdas,lambda_func)
# plot(lambdas,norms)
# abline(v = lambda,lty = 2,col = "red")
# abline(v = -lambda,lty = 2,col ="grey")
