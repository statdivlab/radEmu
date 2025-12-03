# Make lists of constraint functions and their gradients when using a reference taxon

Make lists of constraint functions and their gradients when using a
reference taxon

## Usage

``` r
make_reference_constraints(p, j)
```

## Arguments

- p:

  The number of columns in the design matrix `X`. If you don't know the
  number of columns, you can find it with
  `ncol(radEmu::make_design_matrix(formula = your_formula, data = your_data))`.
  `your_formula` should be the expression you give to `emuFit`'s
  `formula` argument.

- j:

  A single value or a vector of length `p - 1` where `p` is the number
  of columns in the design matrix `X`. If a single value, `j` will be
  used as the reference category for all covariates. If a vector of
  values, `j[k]` will be used as the reference category for the
  covariate in design matrix column `k + 1`.

## Value

A list with elements `constraints_list` and `constraints_grad_list`. The
`constraints_list` is a list of constraint functions for each column `p`
of the design matrix. By default, the constraint for the intercept is
the pseudo Huber median. The constraints for covariates are determined
by reference categories given by the argument `j`. The
`constraints_grad_list` is a list of gradients of each constraint
function.

## Examples

``` r
# two columns in design matrix, reference taxon is taxon 5
list1 <- make_reference_constraints(p = 2, j = 5)

# four columns in design matrix, reference taxon for all covariates is taxon 2
list2 <- make_reference_constraints(p = 4, j = 2)

# four columns in design matrix, reference taxon for covariates 1 and 2 is taxon 3 and
# reference taxon for covariate 3 is taxon 4
list3 <- make_reference_constraints(p = 4, j = c(3, 3, 4))
```
