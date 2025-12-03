# Fit radEmu model. Called by emuFit

Fit radEmu model. Called by emuFit

## Usage

``` r
emuFit_micro(
  X,
  Y,
  B = NULL,
  constraint_fn = NULL,
  maxit = 250,
  tolerance = 1e-05,
  verbose = TRUE,
  warm_start = TRUE,
  c1 = 1e-04,
  max_stepsize = 0.5,
  max_abs_B = 50,
  use_working_constraint = TRUE,
  j_ref = NULL,
  optimize_rows = TRUE,
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

- constraint_fn:

  function g defining constraint on rows of B; g(B_k) = 0 for rows k =
  1, ..., p of B.

- maxit:

  maximum number of coordinate descent cycles to perform before exiting
  optimization

- tolerance:

  tolerance on improvement in log likelihood at which to exit
  optimization

- verbose:

  logical: print information about optimization progress? Default is
  TRUE.

- warm_start:

  logical: begin from "warm start" obtained from linear regression on
  transformed counts? Default is TRUE.

- c1:

  numeric: value of constant in Armijo condition checked in backtracking
  line search

- max_stepsize:

  numeric: maximum sup-norm value of proposed step. Default is 0.5.

- max_abs_B:

  numeric: maximum value elements of B are allowed to take in absolute
  value. Helps prevent optimization failure in larger problems. Defaults
  to 50.

- use_working_constraint:

  logical: set a column of B equal to zero within optimization. Default
  is TRUE.

- j_ref:

  If use_working_constraint is TRUE, column index of column of B to set
  to zero. Default is NULL, in which case this column is chosen to
  maximize the number of nonzero entries of Y_j_ref.

- optimize_rows:

  If use_working_constraint is TRUE, update overall location of rows of
  B relative to column constrained to equal zero under working
  constraint before iterating through updates to columns of B
  individually. Default is TRUE.

- use_discrete:

  If discrete design matrix, use fast discrete implementation.

## Value

A p x J matrix containing regression coefficients (under constraint
g(B_k) = 0)
