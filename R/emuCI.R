#' Construct confidence intervals for parameters estimated in fitted model 
#'
#' This function takes in a radEmu model object and constructs confidence 
#' intervals. Intervals can be constructed using a variety of robust 
#' standard error estimation methods. 
#'
#' @param emuMod An emuMod object fit with \code{emuFit}.
#' @param type The type of standard error estimation approach. Options are 
#' \code{bayesian_subsample}, \code{sandwich}, and \code{bootstrap_percentile}.
#
#' @return A confidence interval for each \code{beta} parameter estimated
#' in the \code{emuMod}.
#' @author David Clausen
#'
#' @export
emuCI <- function(emuMod,
                  nboot,
                  conf_level = 0.95,
                  rows_of_interest = NULL,
                  parallel = FALSE,
                  blocks = NULL,
                  verbose = FALSE,
                  type = "bootstrap_percentile",
                  nsim = 1000,
                  ncore = 1){
  
  J <- ncol(emuMod$Y)
  
  if(is.null(rows_of_interest)){
    rows_of_interest <- 2:nrow(emuMod$B)
    if(nrow(emuMod$B) ==1){
      stop("When nrow(B) = 1, specify that inference on elements of this row is
            desired. Are you sure this is what you want? (Intercept row frequently has no
            useful interpretation!)")
    }
  }
  
  # bootstrap_t is commented out because it does not provide a ci, only sds
  # if(type == "bootstrap_t"){
  #   boot_sds <- emuSD(emuMod = emuMod,
  #                     nboot = nboot,
  #                     rows_of_interest = rows_of_interest,
  #                     parallel = parallel,
  #                     blocks = blocks,
  #                     verbose = verbose,
  #                     ncore = ncore)
  # }
  
  
  if(type == "bayesian_subsample"){
    cis <- emuBayesSub(emuMod,
                       nboot = nboot,
    )
    return(cis)
  }
  
  
  if(type == "sandwich"){
    message("Computing sandwich-based covariance matrix.")
    sandwich_output <- emuSandwich(emuMod,
                                   nsim = nsim)
    
    uppers <- sandwich_output$uppers
    lowers <- sandwich_output$lowers
    pvals <- sandwich_output$pvals
    
    alpha_level <- 1 - conf_level
    sd_row_index <- 1
    cis <- vector("list",length(rows_of_interest))
    message("Constructing confidence intervals")
    for(B_row in rows_of_interest){
      cis[[sd_row_index]] <-
        data.frame(row = B_row,
                   outcome_index = 1:J,
                   estimate = emuMod$B[B_row,],
                   lower = lowers[(B_row - 1)*J + 1:J],
                   upper = uppers[(B_row - 1)*J + 1:J],
                   pval = pvals[(B_row - 1)*J + 1:J],
                   conf_level = conf_level)
      
      sd_row_index <- sd_row_index + 1
    }
  }
  if(type == "bootstrap_percentile"){
    boot_draws <- emuBoot(emuMod = emuMod,
                          nboot = nboot,
                          parallel = FALSE,
                          blocks = blocks,
                          seed = seed,
                          verbose = verbose,
                          ncore = ncore)
    
    alpha_level <- 1 - conf_level
    sd_row_index <- 1
    cis <- vector("list",length(rows_of_interest))
    message("Constructing confidence intervals")
    for(B_row in rows_of_interest){
      cis[[sd_row_index]] <-
        data.frame(row = B_row,
                   outcome_index = 1:J,
                   estimate = emuMod$B[B_row,],
                   lower = sapply(1:J,
                                  function(j)
                                    quantile(sapply(1:nboot,
                                                    function(d)
                                                      boot_draws[[d]][B_row,j]),
                                             alpha_level/2)),
                   upper = sapply(1:J,
                                  function(j)
                                    quantile(sapply(1:nboot,
                                                    function(d)
                                                      boot_draws[[d]][B_row,j]),
                                             1-alpha_level/2)),
                   pval = NA,
                   conf_level = conf_level)
      
      sd_row_index <- sd_row_index + 1
    }
    
  }
  
  return(do.call(rbind,cis))
  
}
