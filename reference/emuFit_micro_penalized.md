# Fit radEmu model with Firth penalty

Fit radEmu model with Firth penalty

## Usage

``` r
emuFit_micro_penalized(
  X,
  Y,
  B = NULL,
  X_cup = NULL,
  constraint_fn = NULL,
  maxit = 500,
  ml_maxit = 5,
  tolerance = 0.001,
  max_step = 5,
  verbose = TRUE,
  max_abs_B = 250,
  j_ref = NULL,
  use_discrete = TRUE
)
```

## Arguments

- X:

  a p x J design matrix

- Y:

  an n x p matrix of nonnegative observations

- B:

  starting value of coefficient matrix (p x J)

- X_cup:

  design matrix for Y in long format. Defaults to NULL, in which case
  matrix is computed from X.

- constraint_fn:

  function g defining constraint on rows of B; g(B_k) = 0 for rows k =
  1, ..., p of B.

- maxit:

  maximum number of coordinate descent cycles to perform before exiting
  optimization

- ml_maxit:

  numeric: maximum number of coordinate descent cycles to perform inside
  of maximum likelihood fits. Defaults to 5.

- tolerance:

  tolerance on improvement in log likelihood at which to exit
  optimization

- max_step:

  numeric: maximum sup-norm for proposed update steps

- verbose:

  logical: report information about progress of optimization? Default is
  TRUE.

- max_abs_B:

  numeric: maximum allowed value for elements of B (in absolute value).
  In most cases this is not needed as Firth penalty will prevent
  infinite estimates under separation. However, such a threshold may be
  helpful in very poorly conditioned problems (e.g., with many nearly
  collinear regressors). Default is 50.

- j_ref:

  which column of B to set to zero as a convenience identifiability
  during optimization. Default is NULL, in which case this column is
  chosen based on characteristics of Y (i.e., j_ref chosen to maximize
  number of entries of Y_j_ref greater than zero).

- use_discrete:

  If discrete design matrix, use fast discrete implementation.

## Value

A p x J matrix containing regression coefficients (under constraint
g(B_k) = 0)
