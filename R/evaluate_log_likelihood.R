

evaluate_log_likelihood <- function(Y_tilde,
                                    D_tilde,
                                    theta,
                                    weights = NULL){

  if(is.null(weights)){
    weights <- Matrix(1,ncol = nrow(Y_tilde),nrow = 1)
  }
  return(
    as.numeric(
      Matrix(weights,nrow = 1)%*%(Y_tilde*(D_tilde%*%theta) - exp(D_tilde%*%theta))
      )
  )
}
