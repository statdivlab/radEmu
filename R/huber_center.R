huber_center <- function(x,
                         d = 1,
                         limit = 20,
                         newton = TRUE,
                         tolerance = 1e-12){

  if(!newton)
    {return(optim(median(x),
               function(y) huber_loss(x - y,d),
               lower = -limit,
               upper = limit,
               method = "Brent",
               control = list(reltol = 1e-12))$par)
  } else{
      c <- median(x)


      dhuber_dc <- Inf
      while(sum(abs(dhuber_dc))>tolerance){
        outside <- abs(x - c)>d
      dhuber_dc <- -d*sum(sign(x[outside] - c)) -
        sum(x[!outside] - c)
      if(sum(!outside)>0){
      dsqhuber_dcsq <- sum(!outside)
      } else{
        dsqhuber_dcsq <- 1
      }

      c <- c - dhuber_dc/dsqhuber_dcsq

      # print(c)
      # print(sum(abs(dhuber_dc)))
      }

      return(c)


    }
}
