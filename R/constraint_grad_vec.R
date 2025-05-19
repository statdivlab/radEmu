

constraint_grad_vec <- function(constraint_grad_fn,
                                js_used,
                                Bk_constr,
                                j_constr,
                                p){
  
  cg <- constraint_grad_fn(Bk_constr[-j_constr])
  cg_restr <- sapply(js_used,
                          function(ind) cg[ind - as.numeric(ind>j_constr)])

  
  return(cg_restr)
}