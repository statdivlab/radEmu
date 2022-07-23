
emuBoot <- function(emuMod,
                    nboot,
                    parallel = FALSE,
                    blocks = NULL,
                    seed = 0,
                    verbose = FALSE,
                    ncore = NULL){


 if(is.null(blocks)){
   blocks <- 1:nrow(emuMod$Y)
 } else{
   if(length(blocks) != nrow(emuMod$Y)){
     stop("Length of blocks must equal number of rows of Y.")
   }
 }
  block_row <- data.frame("block" = blocks,
                          "row" = 1:nrow(emuMod$Y))
  uni_blocks <- unique(blocks)

  if(!parallel){
    boot_results <- vector("list",nboot)
  for(boot_iter in 1:nboot){
    boot_blocks <- sample(uni_blocks,length(uni_blocks),replace = TRUE)
    boot_rows <- do.call(c,lapply(boot_blocks,
                        function(x)
                          with(block_row,
                               row[block ==x])))

    boot_results[[boot_iter]] <-
      with(emuMod,
           emuFit(X = X[boot_rows,],
             Y = Y[boot_rows,],
             B = B,
             tolerance = tolerance,
             maxit = maxit,
             verbose = verbose,
             constraint_fn = constraint_fn,
             maxit_glm = maxit_glm,
             method = method,
             reweight = FALSE,
             weights = weights[boot_rows,],
             return_a_lot = FALSE)$B)
  }
}

  if(parallel){

    if(!is.null(ncore)){
    boot_results <-
      parallel::mclapply(1:nboot,
                         function(x){
                           set.seed(x + seed);
    boot_blocks <- sample(uni_blocks,length(uni_blocks),replace = TRUE);
    boot_rows <- do.call(c,lapply(boot_blocks,
                                  function(x)
                                    with(block_row,
                                         row[block ==x])));
      return(with(emuMod,
           emuFit(X = X[boot_rows,],
                  Y = Y[boot_rows,],
                  B = B,
                  tolerance = tolerance,
                  maxit = maxit,
                  verbose = verbose,
                  constraint_fn = constraint_fn,
                  maxit_glm = maxit_glm,
                  method = method,
                  reweight = FALSE,
                  weights = weights[boot_rows,],
                  return_a_lot = FALSE)$B))}
      )} else{
        boot_results <-
          parallel::mclapply(1:nboot,
                             function(x){
                               set.seed(x + seed);
                               boot_blocks <- sample(uni_blocks,length(uni_blocks),replace = TRUE);
                               boot_rows <- do.call(c,lapply(boot_blocks,
                                                             function(x)
                                                               with(block_row,
                                                                    row[block ==x])));
                               return(with(emuMod,
                                           emuFit(X = X[boot_rows,],
                                                  Y = Y[boot_rows,],
                                                  B = B,
                                                  tolerance = tolerance,
                                                  maxit = maxit,
                                                  verbose = verbose,
                                                  constraint_fn = constraint_fn,
                                                  maxit_glm = maxit_glm,
                                                  method = method,
                                                  reweight = FALSE,
                                                  weights = weights,
                                                  return_a_lot = FALSE)$B))},
                             mc.cores = ncore)
      }
  }

  return(boot_results)


}
