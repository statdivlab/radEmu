

hess_pseudohuber_center <- function(x,
                                    d,
                                    ind_1,
                                    ind_2){

  ps_center <- pseudohuber_center(x,d)
  scaled_sq <- ((x - ps_center)/d)^2
  w <- sqrt(1/(1 + scaled_sq))
  w3 <- w^3
  w3_sum <- sum(w3)

  w_derivs_d_ind_2 <- -w3*((x - ps_center)/(d^2))*(as.numeric((1:length(x)) == ind_2) - w3[ind_2]/w3_sum)


  return((w3[ind_1]/w3_sum)*(
    (3/w[ind_1])*w_derivs_d_ind_2[ind_1] -
      (1/w3_sum)*sum((3*w^2)*w_derivs_d_ind_2)
  ))

}
