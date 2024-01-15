test_that("in saturated case, augmentations reduce to haldane correction", {
    Y <- matrix(1:4,nrow = 2)
    X <- cbind(1,c(0,1))
    J <- ncol(Y)
    n <- nrow(Y)


    ml_fit <- emuFit_micro(X,Y ,
                           constraint_fn = function(x) x[2],
                           tolerance = 1e-5,
                           verbose = FALSE)

    X_cup <- X_cup_from_X(X,J)
    G <- get_G_for_augmentations(X = X,J = J,n = n,
                                 X_cup = X_cup)
    expect_true(max(abs(get_augmentations(X = X,
                                          G = G,
                                          Y = Y,
                                          B = ml_fit) - 0.5)) < 1e-10)

    #correction should not depend on B in this case -- check this too
    expect_true(max(abs(get_augmentations(X = X,
                                          G = G,
                                          Y = Y,
                                          B = 0*ml_fit) - 0.5)) < 1e-10)

  })


# test_that("with continuous predictor, augmentations match what we get from
#           multinomial directly", {
#             skip("Requires brglm2")
#   set.seed(4323)
#   X <- cbind(1,rnorm(10))
#   z <- rnorm(10) +5
#   J <- 10
#   p <- 2
#   n <- 10
#   b0 <- rnorm(J)
#   b1 <- seq(1,5,length.out = J)
#   b1 <- b1 - mean(b1)
#   b <- rbind(b0,b1)
#   Y <- matrix(NA,ncol = J, nrow = n)
# 
#   for(i in 1:n){
#     for(j in 1:J){
#       temp_mean <- exp(X[i,,drop = FALSE]%*%b[,j,drop = FALSE] + z[i])
#       Y[i,j] <- rnbinom(1, mu= temp_mean,size = 2)*rbinom(1,1,0.4)
#       # Y[i,j] <- rpois(1,lambda = temp_mean)
#     }
#   }
# 
#   pl_fit_new <- emuFit_micro_penalized(X,
#                                        Y,
#                                        B = NULL,
#                                        # constraint_fn = function(x) mean(x),
#                                        maxit = 1000,
#                                        tolerance = 1e-5,
#                                        verbose= TRUE)
# 
#   # sum(pl_fit_new$Y_augmented - Y)
# 
#   x <- X[,2]
# 
#   # Construct a variable with the multinomial categories according to
#   # the HCV and nonABC columns
#   library(brglm2)
#   # data("hepatitis", package = "brglm2")
#   #
#   # hepat <- hepatitis
#   # hepat$type <- with(hepat, factor(1 - HCV * nonABC + HCV + 2 * nonABC))
#   # hepat$type <- factor(hepat$type, labels = c("noDisease", "C", "nonABC"))
#   # contrasts(hepat$type) <- contr.treatment(3, base = 1)
#   #
#   # # Maximum likelihood estimation fails to converge because some estimates are infinite
#   # hepML <- brmultinom(type ~ group * time, data = hepat, weights = counts, type = "ML", slowit = 0.1)
#   # hep_meanBR <- brmultinom(type ~ group * time, data = hepat, weights = counts, type = "AS_mean")
#   # brglm2::brmultinom(Y~x,
#   #                    data = data.frame(x = x))
#   #
#   #
#   colnames(Y) <- paste("tax",1:J,sep = "_")
#   Y_fac <- factor(rep(colnames(Y),each = 10),
#                   levels = paste("tax",1:J,sep = "_"))
#   weights <- do.call(c,lapply(1:J,function(i) Y[,i]))
#   longer_dat <- data.frame(Y = Y_fac,
#                      weights = weights,
#                      x = x)
#   meanBR <- brmultinom(Y ~ x, data = longer_dat, weights = weights,
#                            type = "AS_mean",
#                        control = list(trace = TRUE))
# 
#   plot(  pl_fit_new$B[2,],
#          c(0,summary(meanBR)$coefficients[,2]))
# 
#   brglm_coefs <- c(0,summary(meanBR)$coefficients[,2])
#   names(brglm_coefs) <- NULL
# 
#   pl_coefs <- pl_fit_new$B[2,] - pl_fit_new$B[2,1]
# 
#   expect_equal(pl_coefs,
#               brglm_coefs,
#               tolerance = 1e-4
#      )
# 
# 
# 
# })
