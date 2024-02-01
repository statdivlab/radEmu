

#derivative of pseudohuber center wrt its arguments
#used for testing
#i.e., this gives us little h described in inference 
#section of Clausen & Willis (2024)
dpseudohuber_center_dx <- function(x,
                                   d = 1){
  ps_center <- pseudohuber_center(x,d)
  scaled_sq <- ((x - ps_center)/d)^2
  #derivation of why the following is the derivative of 
  #the pseudohuber centering is given in supplement of 
  #Willis & Clausen (2024)
  w <- sqrt(1/(1 + scaled_sq))
  w3 <- w^3
  return((w3/sum(w3)))
}
