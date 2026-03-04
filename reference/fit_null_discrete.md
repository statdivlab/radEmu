# fits model with B_kj constrained to equal g(B_k) for constraint fn g, for a symmetric constraint with a discrete design

fits model with B_kj constrained to equal g(B_k) for constraint fn g,
for a symmetric constraint with a discrete design

## Usage

``` r
fit_null_discrete(
  Y,
  X,
  k_constr,
  j_constr,
  j_ref,
  trackB = FALSE,
  maxit = 5000,
  tol = 0.1,
  verbose = FALSE,
  ls_max = 20,
  ls_rho = 0.5,
  max_step = 1,
  constraint = "pseudohuber"
)
```

## Arguments

- Y:

  Y (with augmentations)

- X:

  design matrix

- k_constr:

  row index of B to constrain

- j_constr:

  col index of B to constrain

- j_ref:

  column index of convenience constraint

- trackB:

  track value of beta across iterations and return? Default is `FALSE`.

- maxit:

  maximum iterations. Default is `5000`.

- tol:

  tolerance for stopping criteria. Algorithm stops when the root mean of
  the norm of the score vector is less than the tolerance. Default is
  `0.01`.

- verbose:

  should the algorithm print updates for you? Default is `FALSE`.

- ls_max:

  maximum number of iterations in the line search. Default is `20`.

- ls_rho:

  scaling factor in the line search. Default is `0.5`.

- max_step:

  step capping after the line sesarch. Default is `1`.

- constraint:

  What type of symmetric constraint do we have? Options are `"mean"` and
  `"pseudohuber`.

## Value

A list containing elements `B`, `k_constr`, `j_constr`, `niter` and
`Bs`. `B` is a matrix containing parameter estimates under the null
(obtained by maximum likelihood on augmented observations Y),
`k_constr`, and `j_constr` give row and column indexes of the parameter
fixed to be equal to the constraint function \\g()\\ under the null.
`niter` is a scalar giving total number of outer iterations used to fit
the null model, and `Bs` is a data frame containing values of B by
iteration if `trackB` was set equal to TRUE (otherwise it contains a
NULL value).
