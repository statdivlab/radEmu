

emuFit_one <- function(Y,
                       X,
                       rect_weights,
                       j,
                       B,
                       z,
                       method,
                       maxit_glm,
                       WD,
                       info_inv){
  n <- nrow(Y)
  J <- ncol(Y)
  if(method == "ML"){

    if(length(unique(rect_weights[,j])) ==1){
      # glmfit <- suppressWarnings(try(glm(y~X - 1,
      #                                    offset = z,
      #                                    family = "poisson",
      #                                    # start = B[,j],
      #                                    # mustart =
      #                                    # weights = rect_weights[,j],
      #                                    control = list(
      #                                      maxit = maxit_glm),
      #                                    data = data.frame(y = Y[,j]
      #                                    )),
      #                                silent = TRUE))

      glmfit <- suppressWarnings(try(fastglm::fastglm(x = X,
                                                      offset = z,
                                                      family = "poisson",
                                                      # start = B[,j],
                                                      # mustart =
                                                      # weights =
                                                      #   rect_weights[,j]/mean(rect_weights[,j]),
                                                      control = list(
                                                        maxit = maxit_glm),
                                                      y = Y[,j]
      ),
      silent = TRUE))

      glmfit$coef <- glmfit$coefficients
    } else{
    # glmfit <- suppressWarnings(try(glm(y~X - 1,
    #                                    offset = z,
    #                                    family = "poisson",
    #                                    # start = B[,j],
    #                                    # mustart =
    #                                    weights =
    #                                      rect_weights[,j]/mean(rect_weights[,j]),
    #                                    control = list(
    #                                      maxit = maxit_glm),
    #                                    data = data.frame(y = Y[,j]
    #                                    )),
    #                                silent = TRUE))
      glmfit <- suppressWarnings(try(fastglm::fastglm(x = X,
                                         offset = z,
                                         family = "poisson",
                                         # start = B[,j],
                                         # mustart =
                                         weights =
                                           rect_weights[,j]/mean(rect_weights[,j]),
                                         control = list(
                                           maxit = maxit_glm),
                                         y = Y[,j]
                                         ),
                                     silent = TRUE))

      glmfit$coef <- glmfit$coefficients
    }


    if(!is.list(glmfit)){
      glmfit <- suppressWarnings(try(glm(y~X - 1,
                                         offset = z,
                                         family = "poisson",
                                         # start = B[,j],
                                         # weights = rect_weights[,j],
                                         control = list(
                                           maxit = maxit_glm),
                                         data = data.frame(y = Y[,j]
                                         )),
                                     silent = TRUE))

      glmfit <-
        suppressWarnings(try(glm(y~X - 1,
                                 offset = z,
                                 family = "poisson",
                                 start = glmfit$coef,
    weights =
      rect_weights[,j]/mean(rect_weights[,j]),
                                 control = list(
                                   maxit = maxit_glm),
                                 data = data.frame(y = Y[,j]
                                 )),
                             silent = TRUE))
    }




  }
  # if(method == "FL"){
  #   rel_indices <- sapply(1:n, function(i) (i - 1)*J + j)
  #   WD_relevant <- WD[rel_indices,,drop = FALSE]
  #   # Y_for_glm[,j] <- Y[,j]
  #   augmentation <- numeric(n)
  #   augmentation <-
  #     Matrix::diag(Matrix::tcrossprod(Matrix::tcrossprod(WD_relevant,
  #                                                        info_inv),
  #                                     WD_relevant))
  #
  #
  #   Y[,j] <- Y[,j] + augmentation
  #
  #
  #   glmfit <- try(suppressWarnings(glm(y~X - 1,
  #                                      offset = z,
  #                                      start =   B[,j],
  #                                      control = list(maxit = maxit_glm),
  #                                      family = poisson(link = "log"),
  #                                      weights = rect_weights[,j],
  #                                      data = data.frame(y = Y[,j]
  #                                      ))),
  #                 silent= TRUE)
  #   if(!is.list(glmfit)){
  #     glmfit <- suppressWarnings(try(glm(y~X - 1,
  #                                        offset = z,
  #                                        # start =   B[,j],
  #                                        control = list(maxit = maxit_glm),
  #                                        weights = rect_weights[,j],
  #                                        family = poisson(link = "log"),
  #                                        data = data.frame(y = Y[,j]
  #                                        ))))
  #   }
  #
  # }
  if(is.list(glmfit)){
  return(glmfit$coef)
  } else{
    return(NULL)
  }
}
