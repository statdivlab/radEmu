get_sim_bs <- function(J){
  evens <- ((1:J)%%2 ==0)
  b0 <- numeric(J)
  b0[!evens] <- seq(-3,3,length.out = sum(!evens))
  b0[evens] <- seq(3,-3,length.out = sum(evens))
  b1 <- 5*sinh(seq(-10,10,length.out= J))/sinh(10)
  b1[(J/2):(J/2 + 1)] <- 0

  return(list(b0 = b0,
              b1 = b1))
}
