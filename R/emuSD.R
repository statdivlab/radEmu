
emuSD <- function(emuMod,
                  nboot,
                  rows_of_interest = NULL,
                  parallel = FALSE,
                  blocks = NULL,
                  verbose = FALSE,
                  ncore = NULL){

  if(is.null(rows_of_interest)){
    rows_of_interest <- 1:nrow(emuMod$B)
  }

  J <- ncol(emuMod$B)

  boot_results <- emuBoot(emuMod,
                          nboot,
                          parallel = parallel,
                          blocks = blocks,
                          verbose = verbose,
                          ncore = ncore)


  sds <- vector("list",
                length(rows_of_interest))

  row_counter <- 1
  for(k in rows_of_interest){
    sds[[row_counter]] <-
      apply(do.call(rbind,
              lapply(1:nboot, function(x) boot_results[[x]][k,,drop = FALSE])),
            2,sd)

    names(sds[[row_counter]]) <- paste("outcome",1:J,sep = "_")

    names(sds)[row_counter] <- paste("B_row_",k,collapse = "",sep = "")

    row_counter <- row_counter + 1

  }

  return(do.call(rbind,sds))
}
