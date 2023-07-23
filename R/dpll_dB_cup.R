
dpll_dB_cup <- function(X,Y,B){
  J <- ncol(Y)
  dB = B*0
  for(j in 1:J){
    dB[,j] = dd_pll_dBj(X,
                        Y,
                        B,
                        j)
  }
  dB_long = B_cup_from_B(dB)
  return(dB_long)
}
