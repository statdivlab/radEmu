test_that("in saturated case, augmentations reduce to haldane correction", {
    Y <- matrix(1:4,nrow = 2)
    X <- cbind(1,c(0,1))
    J <- ncol(Y)
    n <- nrow(Y)


    ml_fit <- emuFit_micro(X,Y ,
                           constraint_fn = rep(list(function(x) x[2]), 2), 
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


#
