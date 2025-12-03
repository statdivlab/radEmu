# fits model with B_kj constrained to equal g(B_k) for constraint fn g, for a single category constraint

fits model with B_kj constrained to equal g(B_k) for constraint fn g,
for a single category constraint

## Usage

``` r
fit_null_scc(
  B,
  Y,
  X,
  X_cup = NULL,
  k_constr,
  j_constr,
  j_ref,
  constraint_fn,
  constraint_grad_fn,
  rho_init = 1,
  tau = 1.2,
  kappa = 0.8,
  B_tol = 0.01,
  inner_tol = 0.01,
  constraint_tol = 1e-04,
  max_step = 5,
  c1 = 1e-04,
  maxit = 1000,
  inner_maxit = 25,
  verbose = FALSE,
  trackB = FALSE
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

- rho_init:

  where to start quadratic penalty parameter

- tau:

  how much to increment rho by each iteration

- kappa:

  cutoff above which to increment rho. If distance to feasibility
  doesn't shrink by at least this factor in an iteration, increment rho
  by tau.

- B_tol:

  tolerance for convergence in \\max\_{k,j} \lvert B^t\_{kj} - B^{(t -
  1)}\_{kj}\rvert\\

- inner_tol:

  tolerance for inner loop

- constraint_tol:

  tolerance for \\\lvert B_kj - g(B_k)\rvert \\

- max_step:

  maximum step size

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

## Value

A list containing elements `B`, `k_constr`, `j_constr`, `niter` `gap`,
`u`, `rho`, and `Bs`. `B` is a matrix containing parameter estimates
under the null (obtained by maximum likelihood on augmented observations
Y), `k_constr`, and `j_constr` give row and column indexes of the
parameter fixed to be equal to the constraint function \\g()\\ under the
null. `niter` is a scalar giving total number of outer iterations used
to fit the null model, `gap` gives the final value of \\g(B\_{k
constr}) - B\_{k constr, j constr}\\, `u` and `rho` are final values of
augmented Lagrangian parameters, and `Bs` is a data frame containing
values of B by iteration if `trackB` was set equal to TRUE (otherwise it
contains a NULL value). - update based on new algorithm
