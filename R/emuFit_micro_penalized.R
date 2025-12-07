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
#' @param use_discrete If discrete design matrix, use fast discrete implementation.
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
    #get design matrix we'll use for computing augmentations

    if (verbose) {
      message(
        "Constructing expanded design matrix. For larger datasets this
may take a moment."
      )
    }
    if (is.null(X_cup)) {
      X_cup <- X_cup_from_X(X, J)
    }
    G <- get_G_for_augmentations(X, J, n, X_cup)

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
