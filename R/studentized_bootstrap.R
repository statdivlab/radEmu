

studentized_bootstrap <- function(emuMod,
                                  nboot_outer,
                                  nboot_inner,
                                  blocks = NULL,
                                  alpha = 0.05,
                                  seed = 0,
                                  ncore = 5,
                                  in_parallel = FALSE){


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

    boot_results_est <- vector("list",nboot_outer)
    boot_results_se <- vector("list",nboot_outer)
  if(!in_parallel){
      for(boot_iter in 1:nboot_outer){
        print(boot_iter)
        boot_blocks <- sample(uni_blocks,length(uni_blocks),replace = TRUE)
        boot_rows <- do.call(c,lapply(boot_blocks,
                                      function(x)
                                        with(block_row,
                                             row[block ==x])))

        temporary_mod <- emuMod
        temporary_mod$Y <- Y[boot_rows,]
        temporary_mod$X <- X[boot_rows,]
        temporary_mod$weights <- emuMod$weights[boot_rows,]

        boot_results_se[[boot_iter]] <-
          emuSD(emuMod = temporary_mod,
                          nboot = nboot_inner,
                          blocks = boot_blocks)

        boot_results_est[[boot_iter]] <-
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
  } else{
     boot_results <- parallel::mclapply(1:nboot_outer,
                 function(k){
                   # print(boot_iter)
                   set.seed(k + seed)
                   boot_blocks <- sample(uni_blocks,length(uni_blocks),replace = TRUE)
                   boot_rows <- do.call(c,lapply(boot_blocks,
                                                 function(x)
                                                   with(block_row,
                                                        row[block ==x])))

                   temporary_mod <- emuMod
                   temporary_mod$Y <- Y[boot_rows,]
                   temporary_mod$X <- X[boot_rows,]
                   temporary_mod$weights <- emuMod$weights[boot_rows,]

                   se <-
                     emuSD(emuMod = temporary_mod,
                           nboot = nboot_inner,
                           blocks = boot_blocks)

                   est <-
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

                   return(list("se" = se,
                               "est" = est))
                 },
                 mc.cores = ncore)
     for(boot_iter in 1:nboot_outer){
       boot_results_est[[boot_iter]] <- boot_results[[boot_iter]]$est
       boot_results_se[[boot_iter]] <- boot_results[[boot_iter]]$se
     }
  }



      boot_ts <-
        lapply(1:nboot_outer,
             function(k)
               (boot_results_est[[k]] - emuMod$B)/(boot_results_se[[k]]))

      upper_t <- 0*emuMod$B
      lower_t <- 0*emuMod$B
      outer_boot_ses <- 0*emuMod$B

      boot_cis <- data.frame(row = numeric(0),
                             outcome_index = numeric(0),
                             estimate = numeric(0),
                             lower = numeric(0),
                             upper = numeric(0),
                             pval = numeric(0))

      for(k in 1:nrow(upper_t)){
        for(j in 1:ncol(upper_t)){
          Zs <- sapply(1:nboot_outer,
                       function(d) boot_ts[[d]][k,j])
          upper_t[k,j] <- quantile(Zs,
                                   1 - alpha/2)
          lower_t[k,j] <- quantile(Zs,
                                   alpha/2)
          outer_boot_ses[k,j] <- sd(sapply(1:nboot_outer,
                                           function(d) boot_results_est[[d]][k,j]))
          tstat <- emuMod$B[k,j]/outer_boot_ses[k,j]
          boot_cis <- rbind(boot_cis,
                            data.frame(row = k,
                                       outcome_index = j,
                                       estimate = emuMod$B[k,j],
                                       lower = emuMod$B[k,j] -
                                         outer_boot_ses[k,j]*upper_t[k,j],
                                       upper = emuMod$B[k,j] -
                                         outer_boot_ses[k,j]*lower_t[k,j],
                                       pval = mean(abs(Zs)>=abs(tstat))))
        }
      }



      boot_cis$pval[boot_cis$pval ==0] <- paste("<",signif(1/nboot_outer,2),
                                                collapse = "",
                                                sep = "")

      return(boot_cis)

}
