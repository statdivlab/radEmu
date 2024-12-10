
l_aug <- function(Bj,
                  Bj_constr_no_k_constr,
                  z,
                  Y,
                  X,
                  Bk_constr_no_j_j_constr,
                  k_constr,
                  j,
                  j_constr,
                  constraint_fn){
  #update parameter determined by null constraint
  Bk_constr_j_constr <- 
    constraint_fn(c(Bk_constr_no_j_j_constr,
                    Bj[k_constr]))
  #create object to store updated two columns of B (j and j_constr) in
  Bj_j_constr <- cbind(Bj,0)
  
  #update object with elements of B_j_constr not subject to null constraint 
  Bj_j_constr[-k_constr,2] <- Bj_constr_no_k_constr
  
  #update object with element of B_j_constr  subject to null constraint 
  Bj_j_constr[k_constr,2] <- Bk_constr_j_constr 
  
  #construct log means for Yj and Yj_constr
  log_means <- X%*%Bj_j_constr + matrix(z,ncol = 1)%*%matrix(1,nrow = 1, ncol = 2)
  
  #compute ll
  lil_ll <- sum(Y[,c(j,j_constr)]*log_means - exp(log_means))
  
  #return ll
  return(lil_ll)
  
}
