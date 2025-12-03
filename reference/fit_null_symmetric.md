# fits model with B_kj constrained to equal g(B_k) for constraint fn g, for a symmetric constraint

fits model with B_kj constrained to equal g(B_k) for constraint fn g,
for a symmetric constraint

## Usage

``` r
fit_null_symmetric(
  B,
  Y,
  X,
  X_cup = NULL,
  k_constr,
  j_constr,
  j_ref,
  constraint_fn,
  constraint_grad_fn,
  B_tol = 0.01,
  inner_tol = 0.01,
  c1 = 0.01,
  maxit = 1000,
  inner_maxit = 25,
  verbose = FALSE,
  trackB = FALSE,
  use_optim = FALSE,
  ignore_stop = FALSE,
  tol_lik = 1e-05,
  tol_test_stat = 0.01,
  null_window = 5,
  max_step = 1
)
```

## Arguments

- B:

  description

- Y:

  Y (with augmentations)

- X:

  design matrix

- X_cup:

  design matrix for Y in long format. Defaults to NULL, in which case
  matrix is computed from X.

- k_constr:

  row index of B to constrain

- j_constr:

  col index of B to constrain

- j_ref:

  column index of convenience constraint

- constraint_fn:

  constraint function

- constraint_grad_fn:

  gradient of constraint fn

- B_tol:

  tolerance for convergence in \\max\_{k,j} \lvert B^t\_{kj} - B^{(t -
  1)}\_{kj}\rvert\\

- inner_tol:

  tolerance for inner loop

- c1:

  constant for armijo rule

- maxit:

  maximum iterations

- inner_maxit:

  max iterations per inner loop

- verbose:

  shout at you?

- trackB:

  track value of beta across iterations and return?

- use_optim:

  whether to use `optim` instead of fisher scoring. Default is FALSE.

- ignore_stop:

  whether to ignore stopping criteria and run `maxit` iterations (could
  be helpful for diagnostic plots).

- tol_lik:

  tolerance for relative changes in likelihood for stopping criteria.
  Default is `1e-5`.

- tol_test_stat:

  tolerance for relative changes in test statistic for stopping
  criteria. Default is `0.01`.

- null_window:

  window to use for stopping criteria (this many iterations where
  stopping criteria is met). Default is `5`.

- max_step:

  Default is `1`

## Value

A list containing elements `B`, `k_constr`, `j_constr`, `niter` and
`Bs`. `B` is a matrix containing parameter estimates under the null
(obtained by maximum likelihood on augmented observations Y),
`k_constr`, and `j_constr` give row and column indexes of the parameter
fixed to be equal to the constraint function \\g()\\ under the null.
`niter` is a scalar giving total number of outer iterations used to fit
the null model, and `Bs` is a data frame containing values of B by
iteration if `trackB` was set equal to TRUE (otherwise it contains a
NULL value). - update based on new algorithm
