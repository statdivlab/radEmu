

generate_sim_setting <- function(J,
                                 n_correlations = NULL,
                                 max_no_correlated = NULL

){

  if(is.null(n_correlations)){
    n_correlations <- ceiling(J/3)
  }
  if(is.null(max_no_correlated)){
    max_no_correlated <- ceiling(J/10)
  }
  if(max_no_correlated ==1){
    stop("Maximum number of taxa in each correlation term must be >1")
  }

  cor_structs <- lapply(1:n_correlations,
                        function(x) sample(1:J,sample(2:max_no_correlated, size = 1)) %>% # size arg added by Amy 
                          (function(y)
                            as.data.frame(cbind("j" = y,
                                               "sgn" = 1- 2*rbinom(length(y),
                                                                   1,0.5)))))

  return(cor_structs)


}
