yy <- Y
xx <- model.matrix(~ treatment)
xx
yy
#efm <- emuFit_micro(X = xx, Y = yy)
efm

yy
xx

emuFit_micro_discrete <- function(
    xx,
    yy
) {
  
  j_ref <- 1
  n <- nrow(yy)
  J <- ncol(yy)
  p <- ncol(xx)
  
  distinct_xx <- unique(xx)
 
  stopifnot(ncol(xx) == nrow(distinct_xx))
  

  # fitted values eta = X %*% beta
  etahats <- matrix(NA, nrow = p, ncol = J) # p x J
  etahats[, j_ref] <- 0
  
  groups <- split(
    seq_len(nrow(xx)),                        # row indices of original data
    apply(xx, 1, function(r)
      paste(r, collapse = "_"))               # grouping key by row contents
  )
  
  for (the_cat in 1:length(groups)) {
    the_xs <- groups[[the_cat]]
    totals <- apply(yy[the_xs, ], 2, sum)
    pihats <- totals / sum(totals)
    etahats[the_cat, setdiff(1:J, j_ref)] <- log(pihats[-j_ref] / pihats[j_ref])
  }
  # beta = X^{-1} %*% eta
  betahats <- round(ginv(distinct_xx), 8) %*% etahats
  #list(betahats, etahats)
  betahats
}

distinct_xx


groups
xx
