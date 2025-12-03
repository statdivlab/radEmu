# Run robust score test

Run robust score test

## Usage

``` r
score_test(
  B,
  Y,
  X,
  X_cup,
  k_constr,
  j_constr,
  constraint_fn,
  constraint_grad_fn,
  rho_init = 1,
  tau = 2,
  kappa = 0.8,
  B_tol = 0.001,
  inner_tol = 0.01,
  constraint_tol = 0.001,
  j_ref,
  c1 = 1e-04,
  maxit = 1000,
  inner_maxit = 25,
  ntries = 4,
  verbose = FALSE,
  trackB = FALSE,
  I_inv = NULL,
  Dy = NULL,
  return_both_score_pvals = FALSE,
  cluster = NULL,
  null_diagnostic_plots = FALSE,
  ignore_stop = FALSE,
  tol_lik = 1e-05,
  tol_test_stat = 0.01,
  null_window = 5
)
```

## Arguments

- B:

  value of coefficient matrix (p x J) returned by full model fit or
  value of coefficient matrix to start null estimation at given as input
  to emuFit

- Y:

  an n x J matrix or dataframe of *augmented* nonnegative observations
  (i.e., observations Y plus augmentations from last iteration of
  maximum penalized likelihood estimation for full model)

- X:

  an n x p matrix of covariates

- X_cup:

  the design matrix for long format Y in long format B (nJ x pJ)

- k_constr:

  row index of element of B to be tested for equality to row
  identifiability constraint

- j_constr:

  column index of element of B to be tested for equality to row
  identifiability constraint

- constraint_fn:

  function g defining a constraint on rows of B; g(B_k) = 0 for rows k =
  1, ..., p of B. Default function is a smoothed median (minimizer of
  pseudohuber loss).

- constraint_grad_fn:

  function returning gradient of constraint function (as a function of a
  row of B)

- rho_init:

  numeric: value at which to initiate rho parameter in augmented
  Lagrangian algorithm. Default is 1.

- tau:

  numeric: value to scale rho by in each iteration of augmented
  Lagrangian algorithm that does not move estimate toward zero
  sufficiently. Default is 2.

- kappa:

  numeric: value between 0 and 1 that determines the cutoff on the ratio
  of current distance from feasibility over distance in last iteration
  triggering scaling of rho. If this ratio is above kappa, rho is scaled
  by tau to encourage estimate to move toward feasibility.

- B_tol:

  numeric: convergence tolerance for null model fits for score testing
  (if max of absolute difference in B across outer iterations is below
  this threshold, we declare convergence). Default is 0.001.

- inner_tol:

  numeric: convergence tolerance for inner loop of null fitting
  algorithm (if max of absolute difference in B across inner iterations
  is below this threshold, we declare convergence). Default is 0.01.

- constraint_tol:

  numeric: constraint tolerance for fits under null hypotheses (tested
  element of B must be equal to constraint function to within this
  tolerance for a fit to be accepted as a solution to constrained
  optimization problem). Default is 1e-5.

- j_ref:

  column index of convenience constraint

- c1:

  numeric: parameter for Armijo line search. Default is 1e-4.

- maxit:

  maximum number of outer iterations of augmented lagrangian algorithm
  to perform before exiting optimization. Default is 1000.

- inner_maxit:

  maximum number of coordinate descent passes through columns of B to
  make within each outer iteration of augmented lagrangian algorithm
  before exiting inner loop

- ntries:

  numeric: total number of times to attempt optimization under null if
  optimization fails (optimization parameters will be tweaked in
  subsequent fits to attempt to avoid failure). Default is 4.

- verbose:

  provide updates as model is being fitted? Defaults to TRUE.

- trackB:

  store and return values of B at each iteration of optimization
  algorithm? Useful for debugging. Default is FALSE.

- I_inv:

  Optional: matrix containing inverted information matrix computed under
  full model. Default is NULL, in which case information is recomputed
  under null, which we recommend.

- Dy:

  Optional: matrix containing empirical score covariance computed under
  full model. Default is NULL, in which case this quantity is recomputed
  under null, which we recommend.

- return_both_score_pvals:

  logical: should score p-values be returned using both information
  matrix computed from full model fit and from null model fits? Default
  is FALSE. This parameter is used for simulations - in any applied
  analysis, type of p-value to be used should be chosen before
  conducting tests.

- cluster:

  a numeric vector giving cluster membership for each row of Y to be
  used in computing GEE test statistics. Default is NULL, in which case
  rows of Y are treated as independent.

- null_diagnostic_plots:

  logical: should diagnostic plots be made for estimation under the null
  hypothesis? Default is `FALSE`.

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

## Value

A list containing elements `score_stat`, `pval`,
`log_pval`,'niter`, `convergence`, `gap`, `u`, `rho`, `tau`, `inner_maxit`, `null_B`, and `Bs`. `score_stat`gives the value of the robust score statistic for $H_0: B_{k_constr,j_constr} = g(B_{k_constr})$.`pval`and`log_pval`are the p-value (on natural and log scales) corresponding to the score statistic (log_pval may be useful when the p-value is very close to zero). `gap`is the final value of $g(B_{k_constr}) - B_{k_constr, j_constr}$ obtained in optimization under the null.`u`and`rho`are final values of augmented Lagrangian parameters returned by null fitting algorithm.`tau`is the final value of`tau`that is used to update the`rho`values and`inner_maxit`is the final maximum number of iterations for the inner optimization loop in optimization under the null, in which B and z parameter values are maximized for specific`u`and`rho`parameters.`null_B`is the value of B returned but the null fitting algorithm.`Bs`is by default`NULL`; if `trackB
= TRUE`, `Bs\` is a data frame containing values of B by outcome
category, covariate, and iteration.
