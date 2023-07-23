get_dg_dBj <- function(B,
                       j,
                       constraint_fn,
                       k_star,
                       # j_star = NULL,
                       constraint_type = "mean",
                       huber_param = NULL,
                       for_testing = FALSE){
  J <- ncol(B)
  p <- nrow(B)
  dg_dBj <- rep(0,p)
  if(constraint_type == "mean"){
    if(!for_testing){
  dg_dBj[k_star] <- 1/(J - 1)
  } else{
    dg_dBj[k_star] <- 1/J
  }
  }

  if(constraint_type == "huber"){
    if(is.null(huber_param)){
      huber_param <- 1
      warning("Argument huber_param not provided; defaulting to value 1.")
    }
    if(abs(B[k_star,j])<= huber_param){
      if(!for_testing){
        dg_dBj[k_star] <- 1/(sum(abs(B[k_star,]) <= huber_param) - 1)
      } else{
    dg_dBj[k_star] <- 1/(sum(abs(B[k_star,]) <= huber_param))
    }
    }
  }



  # if(constraint_type == "pseudohuber"){
  #   dg_dBj[k_star] <- numDeriv::grad(B[k_star,j],
  #                                    function(b){
  #                                      x <- B[k_star,];
  #                                      x[j] <- b;
  #                                      x <- x[-j_star]
  #
  #                                      return(pseudohuber_center(x,d = d,limit =
  #                                                                  max(abs(B[k_star,]))))})
  # }



  return(dg_dBj)
}
