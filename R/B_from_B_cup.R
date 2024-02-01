
#convert B_cup (long format) to B

B_from_B_cup <- function(B_cup, #long format beta
                         J, #number taxa
                         p #predictors incl. intercept
                         ){
  B <- matrix(0,nrow = p, ncol = J)

  #populate columns of B from B_cup
  for(j in 1:J){
    B[,j] <-  B_cup[(j - 1)*p + 1:p]
  }
  return(B)
}
