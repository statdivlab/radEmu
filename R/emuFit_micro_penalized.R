#' Fit radEmu model with Firth penalty
#'
#' @param X a p x J design matrix
#' @param Y an n x p matrix of nonnegative observations
#' @param B starting value of coefficient matrix (p x J)
#' @param X_cup design matrix for Y in long format. Defaults to NULL, in
#' which case matrix is computed from X.
#' @param constraint_fn function g defining constraint on rows of B; g(B_k) = 0
#' for rows k = 1, ..., p of B.
#' @param maxit maximum number of coordinate descent cycles to perform before
#' exiting optimization
#' @param ml_maxit numeric: maximum number of coordinate descent cycles to perform inside
#' of maximum likelihood fits. Defaults to 5.
#' @param tolerance tolerance on improvement in log likelihood at which to
#' exit optimization
#' @param max_step numeric: maximum sup-norm for proposed update steps
#' @param verbose logical: report information about progress of optimization? Default is TRUE.
#' @param max_abs_B numeric: maximum allowed value for elements of B (in absolute value). In
#' most cases this is not needed as Firth penalty will prevent infinite estimates
#' under separation. However, such a threshold may be helpful in very poorly conditioned problems (e.g., with many
#' nearly collinear regressors). Default is 50.
#' @param j_ref which column of B to set to zero as a convenience identifiability
#' during optimization. Default is NULL, in which case this column is chosen based
#' on characteristics of Y (i.e., j_ref chosen to maximize number of entries of
#' Y_j_ref greater than zero).
#' @param use_discrete logical: if the design matrix is discrete (it has exactly p distinct
#' rows), compute the penalized estimate and data augmentations in closed form rather than
#' iteratively. Default is TRUE.
#'
#' @return A p x J matrix containing regression coefficients (under constraint
#' g(B_k) = 0)
#'
emuFit_micro_penalized <-
  function(
    X,
    Y,
    B = NULL,
    X_cup = NULL,
    constraint_fn = NULL,
    maxit = 500,
    ml_maxit = 5,
    tolerance = 1e-3,
    max_step = 5,
    verbose = TRUE,
    max_abs_B = 250,
    j_ref = NULL,
    use_discrete = TRUE
  ) {
    J <- ncol(Y)
    p <- ncol(X)
    n <- nrow(Y)
    Y_augmented <- Y
    if (is.null(B)) {
      fitted_model <- NULL
    } else {
      fitted_model <- B
    }
    converged <- FALSE
    counter <- 0

    #for discrete designs, the penalized estimate is available in closed form:
    #the unpenalized discrete estimate from covariate-pattern-by-category totals
    #with 1/2 added to each total. 
    groups <- discrete_groups(X)
    
    if (use_discrete & nrow(groups$distinct_X) == p) {
      
      if (qr(groups$distinct_X)$rank < p) {
        
        stop(
          "Design matrix X inputted for the model is rank-deficient, preventing proper model fitting.
  This might be due to multicollinearity, overparameterization, or redundant factor levels included in covariates.
  Consider removing highly correlated covariates or adjusting factor levels to ensure a full-rank design. \n"
        )
      }
      if (is.null(constraint_fn)) {
        constraint_fn <- rep(list(function(x) pseudohuber_median(x, 0.1)), p)
      }
      totals <- rowsum(Y, groups$group, reorder = TRUE)
      if (anyNA(totals)) { # integer overflow
        totals <- rowsum(1 * Y, groups$group, reorder = TRUE)
      }
      B <- emuFit_micro_discrete(X = groups$distinct_X,
                                 Y = totals + 0.5,
                                 j_ref = j_ref)
      B[B < -max_abs_B] <- -max_abs_B
      B[B > max_abs_B] <- max_abs_B
      for (k in 1:p) {
        B[k, ] <- B[k, ] - constraint_fn[[k]](B[k, ])
      }
      Y_augmented <- Y + get_augmentations_discrete(X = X, Y = Y, B = B,
                                                    groups = groups)
      return(list(
        "Y_augmented" = Y_augmented,
        "B" = B,
        "convergence" = TRUE
      ))
    }

    #get design matrix we'll use for computing augmentations

    if (verbose) {
      message(
        "Constructing expanded design matrix. For larger datasets this
may take a moment."
      )
    }
    if (is.null(X_cup)) {
      X_cup <- X_cup_from_X_fast(X, J)
    }
    G <- get_G_for_augmentations_fast(X, J, n, X_cup)

    while (!converged) {
      # print(counter)

      if (counter == 0 & is.null(B)) {
        Y_augmented <- Y + 1e-3 * mean(Y) #ensures we don't diverge to
        #infinity in first iteration
        #after which point we use
        #data augmentations based on B
        #is there a smarter way to start?
        #probably.
      } else {
        if (verbose) {
          message(
            "Computing data augmentations for Firth penalty. For larger models, this may take some time."
          )
        }

        augmentations <- get_augmentations(
          X = X,
          G = G,
          Y = Y,
          B = fitted_model
        )
        Y_augmented <- Y + augmentations
      }
      if (!is.null(fitted_model)) {
        old_B <- fitted_model
      } else {
        old_B <- Inf
      }
      #fit model by ML to data with augmentations
      fitted_model <- emuFit_micro(
        X,
        Y_augmented,
        B = fitted_model,
        constraint_fn = constraint_fn,
        maxit = ml_maxit,
        warm_start = TRUE,
        max_abs_B = max_abs_B,
        use_working_constraint = TRUE,
        max_stepsize = max_step,
        tolerance = tolerance,
        verbose = verbose,
        j_ref = j_ref,
        use_discrete = use_discrete
      )

      B_diff <- max(abs(fitted_model - old_B)[abs(fitted_model) < max_abs_B])

      if (B_diff < tolerance) {
        converged <- TRUE
        actually_converged <- TRUE
      }

      if (counter > maxit) {
        converged <- TRUE
        actually_converged <- FALSE
      }
      counter <- counter + 1
    }

    return(list(
      "Y_augmented" = Y_augmented,
      "B" = fitted_model,
      "convergence" = actually_converged
    ))
  }
