

emuBayesSub <- function(emuMod,
                        nboot,
                        m = NULL,
                        rows_of_interest = NULL,
                        parallel = FALSE,
                        ncore = 1,
                        seed = 0,
                        conf_level = 0.95){

  alpha_level <- 1 - conf_level
  n <- nrow(emuMod$Y)
  J <- ncol(emuMod$Y)
  if(is.null(m)){
    m <- sqrt(n)
  }
  p <- nrow(emuMod$B)
  if(is.null(rows_of_interest)){
    rows_of_interest <- 2:p
  }
  if(!parallel){
  boot_results <- vector("list",nboot)
  for(boot_iter in 1:nboot){
    boot_weights <- rgamma(n,
                           shape = m/n,
                           scale = 1)
    boot_weights <- boot_weights/sum(boot_weights)
    weights <- matrix(boot_weights,ncol= 1)%*%matrix(1,nrow =1,
                                                     ncol = J)
    weights <- weights*emuMod$weights
    weights <- n*J*weights/sum(weights)

    boot_results[[boot_iter]] <-
      emuFit(X = emuMod$X,
             Y = emuMod$Y,
             B = emuMod$B,
             tolerance = 1e-1,
             optim_only = TRUE,
             maxit = 100,
             verbose = FALSE,
             constraint_fn = emuMod$constraint_fn,
             maxit_glm = emuMod$maxit_glm,
             method = emuMod$method,
             reweight = FALSE,
             weights = weights,
             test_firth = FALSE,
             return_a_lot = TRUE,
             prefit = TRUE)


  }
  } else{
 boot_results <- parallel::mclapply(1:nboot,
                    function(b){
                      set.seed(seed + b)
        boot_weights <- rgamma(n,
                               shape = m/n,
                               scale = 1);
        boot_weights <- boot_weights/sum(boot_weights)
        weights <- matrix(boot_weights,ncol= 1)%*%matrix(1,nrow =1,
                                                         ncol = J);
        weights <- weights*emuMod$weights;
        weights <- n*J*weights/sum(weights);
        return(
          emuFit(X = emuMod$X,
                 Y = emuMod$Y,
                 B = emuMod$B,
                 tolerance = emuMod$tolerance,
                 optim_only = TRUE,
                 maxit = emuMod$maxit,
                 verbose = FALSE,
                 constraint_fn = emuMod$constraint_fn,
                 maxit_glm = emuMod$maxit_glm,
                 method = emuMod$method,
                 reweight = FALSE,
                 weights = weights,
                 test_firth = FALSE,
                 return_a_lot = FALSE,
                 prefit = TRUE))},
        mc.cores = ncore)
  }

  cis <- vector("list",length(rows_of_interest))
  counter <- 1
  for(B_row in rows_of_interest){

  sub_draws <-
    sqrt(m)*do.call(cbind,
            lapply(1:nboot,
           function(d)
             boot_results[[d]]$B[B_row,] - emuMod$B[B_row,]))

  uppers <-
    emuMod$B[B_row,] + (1/sqrt(n))*apply(sub_draws,
                             1,
                             function(x) quantile(x,1 - alpha_level/2))

  lowers <-
    emuMod$B[B_row,] + (1/sqrt(n))*apply(sub_draws,
                                         1,
                                         function(x) quantile(x,
                                                               alpha_level/2)
                                         )

  pvals <- do.call(c,
                   lapply(1:J,
                  function(j)
                    mean(abs((1/sqrt(n))*sub_draws[j,])>=
                           abs(emuMod$B[B_row,j]))))

  cis[[counter]] <- data.frame(row = B_row,
             outcome_index = 1:J,
             estimate = emuMod$B[B_row,],
             lower = lowers,
             upper = uppers,
             pval = pvals,
             conf_level = conf_level)
  counter <- counter + 1
    }

  return(do.call(rbind,cis))
}
