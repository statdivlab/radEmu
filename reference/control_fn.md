# Create list of control options (to pass to `emuFit()`)

Create list of control options (to pass to
[`emuFit()`](https://statdivlab.github.io/radEmu/reference/emuFit.md))

## Usage

``` r
control_fn(
  control = list(),
  max_step = 1,
  ignore_stop = FALSE,
  use_fullmodel_info = FALSE,
  use_fullmodel_cov = FALSE,
  use_both_cov = FALSE,
  inner_maxit = 25,
  inner_tol = 1,
  c1 = 1e-04,
  trackB = FALSE,
  return_nullB = FALSE,
  return_score_components = FALSE,
  return_both_score_pvals = FALSE,
  B_null_tol = 0.001,
  rho_init = 1,
  tau = 2,
  kappa = 0.8,
  constraint_tol = 1e-05,
  ntries = 4
)
```

## Arguments

- control:

  Current control list (optional), will augment it with missing
  arguments

- max_step:

  Maximum stepsize; update directions computed during estimation (under
  the alternative). Will be rescaled if a step in any parameter exceeds
  this value. Defaults to 1.

- ignore_stop:

  whether to ignore stopping criteria and run `maxit` iterations
  (helpful for diagnostic plots to determine convergence).

- use_fullmodel_info:

  Used in estimation under the null hypothesis. Whether to use
  information matrix from estimation under the alternative hypothesis to
  construct the robust score statistic (instead of information matrix
  from estimation under the null hypothesis). Defaults to `FALSE`.

- use_fullmodel_cov:

  Used in estimation under the null hypothesis. Whether to use
  covariance matrix from estimation under the alternative hypothesis to
  construct the robust score statistic (instead of covariance matrix
  from estimation under the null hypothesis). Defaults to `FALSE`.

- use_both_cov:

  Used in estimation under the null hypothesis. Whether to do score test
  twice, once with covariance matrix under the alternative hypothesis
  and once with covariance matrix under the null hypothesis. Defaults to
  `FALSE`.

- inner_maxit:

  Used in estimation under the null hypothesis. Maximum number of
  iterations within each inner loop of estimation under null hypothesis
  algorithm. Default is `25`.

- inner_tol:

  Used in estimation under the null hypothesis. Convergence tolerance
  within each inner loop of estimation under null hypothesis algorithm.
  Default is `1`.

- c1:

  Used in estimation under the null hypothesis. Parameter for Armijo
  line search. Default is `1e-4`.

- trackB:

  Used in estimation under the null hypothesis. When `TRUE` will track
  the value of `B` in each iteration of optimization algorithm. Defaults
  to `FALSE`.

- return_nullB:

  Used in estimation under the null hypothesis. When `TRUE` will return
  the final value of `B` under each null hypothesis tested. Defaults to
  `FALSE`.

- return_score_components:

  Used in estimation under the null hypothesis. When `TRUE` will return
  the components of the robust score test statistic for each null
  hypothesis tested. Defaults to `FALSE`.

- return_both_score_pvals:

  Used in estimation under the null hypothesis, with `use_both_cov`.
  Defaults to `FALSE`.

- B_null_tol:

  Used in estimation under the null hypothesis, for the augmented
  Lagrangian algorithm. numeric: convergence tolerance for null model
  fits for score testing (if max of absolute difference in B across
  outer iterations is below this threshold, we declare convergence).
  Default is `0.001`.

- rho_init:

  Used in estimation under the null hypothesis, for the augmented
  Lagrangian algorithm. Value at which to initiate rho parameter in
  augmented Lagrangian algorithm. Default is `1`.

- tau:

  Used in estimation under the null hypothesis, for the augmented
  Lagrangian algorithm. Value to scale `rho` by in each iteration of
  augmented Lagrangian algorithm that does not move estimate toward zero
  sufficiently. Default is `2`.

- kappa:

  Used in estimation under the null hypothesis, for the augmented
  Lagrangian algorithm. Value between `0` and `1` that determines the
  cutoff on the ratio of current distance from feasibility over distance
  in last iteration triggering scaling of `rho`. If this ratio is above
  `kappa`, `rho` is scaled by `tau` to encourage estimate to move toward
  feasibility.

- constraint_tol:

  Used in estimation under the null hypothesis, for the augmented
  Lagrangian algorithm. Constraint tolerance for fits under null
  hypotheses (tested element of `B` must be equal to constraint function
  to within this tolerance for a fit to be accepted as a solution to
  constrained optimization problem). Default is `1e-5`.

- ntries:

  Used in estimation under the null hypothesis, for the augmented
  Lagrangian algorithm. The number of times to try optimization.
  Successive tries will change `tau` and `inner_maxit` and retry.

## Value

A list containing control options, to have more control over
optimization algorithms used by `radEmu`. This can be passed into
[`emuFit()`](https://statdivlab.github.io/radEmu/reference/emuFit.md).
