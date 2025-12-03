# Runs checks for appropriate arguments before running `emuFit()`

Runs checks for appropriate arguments before running
[`emuFit()`](https://statdivlab.github.io/radEmu/reference/emuFit.md)

## Usage

``` r
emuFit_check(
  Y,
  X = NULL,
  formula = NULL,
  data = NULL,
  assay_name = NULL,
  cluster = NULL,
  B_null_list = NULL,
  test_kj = NULL,
  match_row_names = TRUE,
  verbose = FALSE,
  remove_zero_comparison_pvals = 0.01,
  unobserved_taxon_error = TRUE,
  constraint_fn,
  constraint_grad_fn,
  constraint_param,
  run_score_tests = TRUE,
  null_fit_alg = "constraint_sandwich"
)
```

## Arguments

- Y:

  an n x J matrix or dataframe of nonnegative observations, or a
  phyloseq object containing an otu table and sample data.

- X:

  an n x p matrix or dataframe of covariates (optional)

- formula:

  a one-sided formula specifying the form of the mean model to be fit

- data:

  an n x p data frame containing variables given in `formula`

- assay_name:

  a string containing the desired assay name within a
  `TreeSummarizedExperiment` object. This is only required if Y is a
  `TreeSummarizedExperiment` object, otherwise this argument does
  nothing and can be ignored.

- cluster:

  a vector giving cluster membership for each row of Y to be used in
  computing GEE test statistics. Default is NULL, in which case rows of
  Y are treated as independent.

- B_null_list:

  list of starting values of coefficient matrix (p x J) for null
  estimation. This should either be a list with the same length as
  `test_kj`. If you only want to provide starting values for some tests,
  include the other elements of the list as `NULL`.

- test_kj:

  a data frame whose rows give coordinates (in category j and
  covariate k) of elements of B to construct hypothesis tests for. If
  `test_kj` is not provided, all elements of B save the intercept row
  will be tested.

- match_row_names:

  logical: Make sure rows on covariate data and response data correspond
  to the same sample by comparing row names and subsetting/reordering if
  necessary.

- verbose:

  provide updates as model is being fitted? Defaults to FALSE. If user
  sets verbose = TRUE, then key messages about algorithm progress will
  be displayed. If user sets verbose = "development", then key messages
  and technical messages about convergence will be displayed. Most users
  who want status updates should set verbose = TRUE.

- remove_zero_comparison_pvals:

  Should score p-values be replaced with NA for zero-comparison
  parameters? These parameters occur for categorical covariates with
  three or more levels, and represent parameters that compare a
  covariate level to the reference level for a category in which the
  comparison level and reference level both have 0 counts in all
  samples. These parameters can have misleadingly small p-values and are
  not thought to have scientifically interesting signals. We recommend
  removing them before analyzing data further. If TRUE, all
  zero-comparison parameter p-values will be set to NA. If FALSE no
  zero-comparison parameter p-values will be set to NA. If a value
  between 0 and 1, all zero-comparison p-values below the value will be
  set to NA. Default is `0.01`.

- unobserved_taxon_error:

  logical: should an error be thrown if Y includes taxa that have 0
  counts for all samples? Default is TRUE.

- constraint_fn:

  function g defining a constraint on rows of B; g(B_k) = 0 for rows k =
  1, ..., p of B. Default function is a smoothed median (minimizer of
  pseudohuber loss). If a number is provided a single category
  constraint will be used with the provided category as a reference
  category. This argument can either be a single constraint function to
  be used for all rows of B, or a list of length p of constraints to be
  used for each row of B.

- constraint_grad_fn:

  derivative of constraint_fn with respect to its arguments (i.e.,
  elements of a row of B). If `constraint_fn` is a list of constraint
  functions, then this argument must also be a list.

- constraint_param:

  If pseudohuber centering is used (this is the default), parameter
  controlling relative weighting of elements closer and further from
  center. (Limit as `constraint_param` approaches infinity is the mean;
  as this parameter approaches zero, the minimizer of the pseudo-Huber
  loss approaches the median.)

- run_score_tests:

  logical: perform robust score testing?

- null_fit_alg:

  Which null fitting algorithm to use `"fisher_scoring"` or
  `"augmented_lagrangian"`. Default and recommended approach is
  `"fisher_scoring"`.

## Value

returns objects `Y`, `X`, `cluster`, and `B_null_list`, which may be
modified by tests, and throw any useful errors, warnings, or messages.
