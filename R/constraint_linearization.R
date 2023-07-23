


constraint_linearization <- function(x,
                                     x0,
                                     constraint_fn_at_x0,
                                     constraint_grad_at_x0){

  return(constraint_fn_at_x0 + sum(constraint_grad_at_x0*(x - x0)))
}
