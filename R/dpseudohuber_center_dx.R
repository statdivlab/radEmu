


dpseudohuber_center_dx <- function(x,
                                   d = 1){
  ps_center <- pseudohuber_center(x,d)
  scaled_sq <- ((x - ps_center)/d)^2
  w <- sqrt(1/(1 + scaled_sq))
  w3 <- w^3
  return((w3/sum(w3)))
}
