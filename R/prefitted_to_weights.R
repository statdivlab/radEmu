

prefitted_to_weights <- function(Y,
                                 prefitted,
                                 weight_stabilization = 1){
  n <- nrow(Y)
  J <- ncol(Y)
  square_resids <- (Y - prefitted)^2
  square_resids_long <- as.numeric(Y_to_Y_tilde(square_resids))
  prefitted_long <- as.numeric(Y_to_Y_tilde(prefitted))
  prefitted_order <- order(prefitted_long)
  reweights_long <- prefitted_long


  fitted_var <- numeric(length(prefitted_long))
  fitted_var[prefitted_order] <-
    monotone::monotone((square_resids_long[prefitted_order] +
                          weight_stabilization)/
                         (prefitted_long[prefitted_order] +
                            weight_stabilization),
                       w= 1/(prefitted_long[prefitted_order]^2 +
                               .1*square_resids_long[prefitted_order] +
                                    weight_stabilization))
                       # w= sqrt(prefitted_long[prefitted_order] +
                       #            weight_stabilization)/
                       #   sqrt(square_resids_long[prefitted_order] +
                       #          weight_stabilization))#+
                         # prefitted_long[prefitted_order]^2 +
                         #           weight_stabilization
                         # )))

  # fitted_var[prefitted_order] <-
  #   monotone::monotone(square_resids_long[prefitted_order] -
  #                         prefitted_long[prefitted_order],
  #                      w = 1/prefitted_long[prefitted_order])

  # plot(log(prefitted_long),log(fitted_var))
  # abline(a = 0, b= 1)


  # hist(fitted_var)
#
#




  #,
                       # w= (sum(1/((prefitted_long + 1)^2))/((prefitted_long + 1)^2)[prefitted_order]))

  # plot(log(prefitted_long[prefitted_order]),
  #   log(fitted_var))
  # abline(a = 0, b= 1)

  # if(n*J > 5000){
  #   iso <- isoreg(x = prefitted_long[prefitted_order],
  #                 y = square_resids_long[prefitted_order])
  #   reweights_long[prefitted_order] <-
  #     pmin((prefitted_long[prefitted_order] + 1)/(iso$yf +1),1)
  # } else{
  #   pavafit <- cir::cirPAVA(x = prefitted_long[prefitted_order] + rnorm(length(prefitted_long),
  #                                                                       sd = 1e-6),
  #                           y = square_resids_long[prefitted_order])
  #   reweights_long[prefitted_order] <-
  #     pmin((prefitted_long[prefitted_order] +1)/(pavafit + 1),1)}

  # reweights_long <-
  # (prefitted_long + weight_stabilization)/(fitted_var + weight_stabilization)

  # reweights_long <- (prefitted_long)/(prefitted_long +
  #                                       pmax(fitted_var,0))
  reweights_long <- 1/pmax(fitted_var,1)


  # reweights_long[order(prefitted_long,decreasing = TRUE)] <-
  #   monotone::monotone(reweights_long[order(prefitted_long,decreasing = TRUE)])


#
  # plot(log(prefitted_long)[prefitted_order],log(reweights_long)[prefitted_order],
  # type = "s")


#
#   reweights_long <- prefitted_long/fitted_var
#
  # plot(log(prefitted_long),log(reweights_long))
#
# sd(log(reweights_long))




  # reweights_long <- pmin(reweights_long,1)
  weights <- Y_tilde_to_Y(reweights_long,J)
  weights <- n*J*weights/sum(weights)
  return(weights)
}
