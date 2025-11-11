emuFit_micro_discrete <- function(
    X,
    Y,
    j_ref = NULL
) {
  
  if (is.null(j_ref)) {
    j_ref <- 1
  }
  n <- nrow(Y)
  J <- ncol(Y)
  p <- ncol(X)
  
  distinct_xx <- unique(X)

  # fitted values eta = X %*% beta
  etahats <- matrix(NA, nrow = p, ncol = J) # p x J
  etahats[, j_ref] <- 0
  
  groups <- split(
    seq_len(nrow(X)),                        # row indices of original data
    apply(X, 1, function(r)
      paste(r, collapse = "_"))               # grouping key by row contents
  )
  
  for (the_cat in 1:length(groups)) {
    the_xs <- groups[[the_cat]]
    totals <- apply(Y[the_xs, , drop = FALSE], 2, sum)
    pihats <- totals / sum(totals)
    etahats[the_cat, setdiff(1:J, j_ref)] <- log(pihats[-j_ref] / pihats[j_ref])
  }
  # beta = X^{-1} %*% eta
  betahats <- round(MASS::ginv(distinct_xx), 8) %*% etahats
  #list(betahats, etahats)
  betahats
}

