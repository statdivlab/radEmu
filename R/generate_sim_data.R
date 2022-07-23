
generate_sim_data <- function(n_per_group,
                              b0,
                              b1,
                              delta,
                              z_mean,
                              nb_size,
                              cor_structs,
                              constant_z = FALSE,
                              return_Y_circ = FALSE){

  J <- length(b0)

  if(!constant_z){
  z <- rnorm(2*n_per_group,z_mean)
  } else{
    z <- rep(z_mean,2*n_per_group)
  }

  X <- cbind(1,rep(c(0,1),each = n_per_group))
  B <- rbind(b0 + delta,b1)

  log_means <- X%*%B
  means <- exp(log_means)

  for(k in 1:length(cor_structs)){
    for(i in 1:(2*n_per_group)){
    cor_factor <- 2*rbeta(1,0.5,0.5)
    anti_cor <- 2 - cor_factor
    pos_j <- cor_structs[[k]]$j[cor_structs[[k]]$sgn == 1]
    neg_j <- cor_structs[[k]]$j[cor_structs[[k]]$sgn == -1]
    means[i,pos_j] <- means[i,pos_j]*cor_factor
    means[i,neg_j] <- means[i,neg_j]*anti_cor
    }
  }
  means <- means*exp(matrix(rnorm(2*n_per_group*J,
                                        mean = -.05,
                                        sd = sqrt(.1)),nrow = 2*n_per_group,
                                  ncol = J))

  if(return_Y_circ){
    return(means)
  }

  means <- means*exp(matrix(z,ncol = 1)%*%matrix(1,nrow= 1,ncol = J))

  W <- apply(means,c(1,2),function(x) rnbinom(1,mu = x,size = nb_size))

  return(W)

}
