calculate_firth_penalty <- function(D_tilde,
                                    W,
                                    n_skip){


  info <- Matrix::crossprod(D_tilde,
                            Matrix::crossprod(W,D_tilde))
  penalty <- as.numeric(0.5*Matrix::determinant(info +
                                                  Matrix::Diagonal(
                                                    x = rep(1e-4,
                                                            nrow(info))
                                                  ),
                                                logarithm = TRUE)$modulus) -
    0.5*(nrow(info) - n_skip)*log(1e-4)

  # vals <- eigen(info)$value

  # penalty <- suppressWarnings(0.5*sum(log(vals[1:(length(vals) - n_skip)]),
  #                    na.rm = TRUE))




  return(penalty)
}
