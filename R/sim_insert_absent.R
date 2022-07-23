

sim_insert_absent <- function(ci_object,
                              which_absent,
                              J){
  n_missing <- length(which_absent)
  for(absent_taxon in which_absent[order(which_absent)]){
    ci_object$outcome_index[ci_object$outcome_index>absent_taxon] <-
      ci_object$outcome_index[ci_object$outcome_index>absent_taxon] + 1

    if(absent_taxon == 1){
      piece1 <- NULL
    } else{
      piece1 <- ci_object[1:(absent_taxon - 1),]
    }
    piece2 <- data.frame(row = 2,outcome_index = absent_taxon,
                         estimate = NA,
                         lower = NA,
                         upper = NA,
                         pval = NA,
                         conf_level = NA)

    if(absent_taxon == J){
      piece3 <- NULL
    } else{
      piece3 <- ci_object[(absent_taxon):(J - n_missing),]
    }
    ci_object <- rbind(piece1,piece2,piece3)
    n_missing <- n_missing - 1
  }
  rownames(ci_object) <- NULL
  # rownames(ci_object) <- paste("outcome",1:J,sep = "_")
  return(ci_object)
}
