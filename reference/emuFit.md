# Fit radEmu model

Fit radEmu model

## Usage

``` r
emuFit(
  Y,
  X = NULL,
  formula = NULL,
  data = NULL,
  assay_name = NULL,
  cluster = NULL,
  constraint_fn = pseudohuber_median,
  constraint_grad_fn = dpseudohuber_median_dx,
  constraint_param = 0.1,
  verbose = FALSE,
  match_row_names = TRUE,
  unobserved_taxon_error = TRUE,
  penalize = TRUE,
  B = NULL,
  fitted_model = NULL,
  refit = TRUE,
  tolerance = 1e-04,
  maxit = 1000,
  alpha = 0.05,
  return_wald_p = FALSE,
  compute_cis = TRUE,
  run_score_tests = TRUE,
  test_kj = NULL,
  null_fit_alg = "constraint_sandwich",
  B_null_list = NULL,
  maxit_null = 1000,
  tol_lik = 1e-05,
  tol_test_stat = 0.01,
  null_window = 5,
  null_diagnostic_plots = FALSE,
  remove_zero_comparison_pvals = 0.01,
  control = NULL,
  ...
)
```

## Arguments

- Y:

  an n x J matrix or dataframe of nonnegative observations, or a
  phyloseq object containing an otu table and sample data.

- X:

  an n x p design matrix (either provide `X` or `data` and `formula`)

- formula:

  a one-sided formula specifying the form of the mean model to be fit
  (used with `data`)

- data:

  an n x p data frame containing variables given in `formula`

- assay_name:

  a string containing the desired assay name within a
  `TreeSummarizedExperiment` object. This is only required if Y is a
  `TreeSummarizedExperiment` object, otherwise this argument can be
  ignored.

- cluster:

  a vector giving cluster membership for each row of Y to be used in
  computing GEE test statistics. Default is NULL, in which case rows of
  Y are treated as independent.

- constraint_fn:

  (Optional) User-provided constraint function, if default behavior of
  comparing log fold-difference parameters to smoothed median over all
  categories is not desired. If a number is provided a single category
  constraint will be used with the provided category as a reference
  category. This argument can either be a single constraint function to
  be used for all rows of B, or a list of length p of constraints to be
  used for each row of B.

- constraint_grad_fn:

  (Optional) User-provided derivative of constraint function, if default
  behavior of comparing log fold-difference parameters to smoothed
  median over all categories is not desired. If `constraint_fn` is a
  list of constraint functions, then this argument must also be a list.
  If `constraint_fn` is a single number, or a list that includes a
  single number, then the corresponding `constraint_grad_fn` can be set
  to `NULL`, and will be appropriately set within the function.

- constraint_param:

  (Optional) If the smoothed median is used as a constraint (this is the
  default), parameter controlling relative weighting of elements closer
  and further from center. (Limit as `constraint_param` approaches
  infinity is the mean; as this parameter approaches zero, the minimizer
  of the pseudo-Huber loss approaches the median.) If constraint
  function is not smoothed median (implemented in
  [`radEmu::pseudohuber_median()`](https://statdivlab.github.io/radEmu/reference/pseudohuber_median.md))
  then this argument will be ignored.

- verbose:

  provide updates as model is being fitted? Defaults to `FALSE`. If user
  sets `verbose = TRUE`, then key messages about algorithm progress will
  be displayed. If user sets `verbose = "development"`, then key
  messages and technical messages about convergence will be displayed.
  Most users who want status updates should set `verbose = TRUE`.

- match_row_names:

  logical: If `TRUE`, make sure rows on covariate data and response data
  correspond to the same sample by comparing row names and
  subsetting/reordering if necessary. Default is `TRUE`.

- unobserved_taxon_error:

  logical: should an error be thrown if Y includes taxa that have 0
  counts for all samples? Default is TRUE.

- penalize:

  logical: should Firth penalty be used? Default is `TRUE`. Used in
  estimation.

- B:

  starting value of coefficient matrix (p x J) for estimation. If not
  provided, B will be initiated as a zero matrix. Used in estimation.

- fitted_model:

  a fitted model produced by a separate call to emuFit; to be provided
  if score tests are to be run without refitting the full unrestricted
  model. Default is `NULL`.

- refit:

  logical: if `B` or `fitted_model` is provided, should full model be
  fit (`TRUE`) or should fitting step be skipped (`FALSE`), e.g., if
  score tests are to be run on an already fitted model. Default is
  `TRUE`.

- tolerance:

  tolerance for stopping criterion in estimation; once no element of B
  is updated by more than this value in a single step, we exit
  optimization. Defaults to `1e-3`. Used in estimation.

- maxit:

  maximum number of outer iterations to perform before exiting
  optimization. Default is `1000`. Used in estimation.

- alpha:

  nominal type 1 error level to be used to construct confidence
  intervals. Default is `0.05` (corresponding to 95% confidence
  intervals)

- return_wald_p:

  logical: return p-values from Wald tests? Default is `FALSE`.

- compute_cis:

  logical: compute and return Wald CIs? Default is `TRUE`.

- run_score_tests:

  logical: perform robust score testing? Default is TRUE.

- test_kj:

  a data frame whose rows give coordinates (in category j and
  covariate k) of elements of B to construct hypothesis tests for. If
  you don't know which indices k correspond to the covariate(s) that you
  would like to test, run the function
  [`radEmu::make_design_matrix()`](https://statdivlab.github.io/radEmu/reference/make_design_matrix.md)
  in order to view the design matrix, and identify which column of the
  design matrix corresponds to each covariate in your model. This
  argument is required when running score tests.

- null_fit_alg:

  Which null fitting algorithm to use for score tests:
  `"constraint_sandwich"` or `"augmented_lagrangian"`. Default and
  recommended approach is `"constraint_sandwich"`, unless `J < 20`.

- B_null_list:

  list of starting values of coefficient matrix (p x J) for null
  estimation for score testing. This should either be a list with the
  same length as `test_kj`. If you only want to provide starting values
  for some tests, include the other elements of the list as `NULL`.

- maxit_null:

  maximum number of outer iterations to perform before exiting
  optimization. Default is `1000`. Used in estimation under null
  hypothesis for score tests.

- tol_lik:

  tolerance for relative changes in likelihood for stopping criteria.
  Default is `1e-5`. Used in estimation under null hypothesis for score
  tests with "constraint_sandwich" algorithm.

- tol_test_stat:

  tolerance for relative changes in test statistic for stopping
  criteria. Default is `0.01`. Used in estimation under null hypothesis
  for score tests with "constraint_sandwich" algorithm.

- null_window:

  window to use for stopping criteria (this many iterations where
  stopping criteria is met). Default is `5`. Used in estimation under
  null hypothesis for score tests with "constraint_sandwich" algorithm.

- null_diagnostic_plots:

  logical: should diagnostic plots be made for estimation under the null
  hypothesis? Default is `FALSE`.

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

- control:

  A list of control parameters, to have more control over estimation and
  hypothesis testing. See `control_fn` for details.

- ...:

  Additional arguments. Arguments matching the names of
  [`control_fn()`](https://statdivlab.github.io/radEmu/reference/control_fn.md)
  options are forwarded to that function and override defaults. Unknown
  arguments are ignored with a warning.

## Value

A list containing elements 'coef', 'B', 'penalized', 'Y_augmented',
'z_hat', 'I', 'Dy', and 'score_test_hyperparams' if score tests are run.
Parameter estimates by covariate and outcome category (e.g., taxon for
microbiome data), as well as optionally confidence intervals and
p-values, are contained in 'coef'. Any robust score statistics and score
test p-values are also included in 'coef'. If there are any
zero-comparison parameters in the model, a column 'zero_comparison' is
also included, which is TRUE for any parameters that compare the level
of a categorical covariate to a reference level for a category with only
zero counts for both the comparison level and the reference level. This
check is currently implemented for an arbitrary design matrix generated
using the `formula` and `data` arguments, and for a design matrix with
no more than one categorical covariate if the design matrix `X` is input
directly. 'B' contains parameter estimates in matrix format (rows
indexing covariates and columns indexing outcome category / taxon).
'penalized' is equal to TRUE f Firth penalty is used in estimation
(default) and FALSE otherwise. 'z_hat' returns the nuisance parameters
calculated in Equation 7 of the radEmu manuscript, corresponding to
either 'Y_augmented' or 'Y' if the 'penalized' is equal to TRUE or
FALSE, respectively. I' and 'Dy' contain an information matrix and
empirical score covariance matrix computed under the full model.
'score_test_hyperparams' contains parameters and hyperparameters related
to estimation under the null, including whether or not the algorithm
converged, which can be helpful for debugging.

## Examples

``` r
# data frame example
data(wirbel_sample_small)
data(wirbel_otu_small)
emuRes <- emuFit(formula = ~ Group, data = wirbel_sample_small, Y = wirbel_otu_small,
                 test_kj = data.frame(k = 2, j = 1), tolerance = 0.01) 
 # here we set large tolerances for the example to run quickly, 
 # but we recommend smaller tolerances in practice

# TreeSummarizedExperiment example (only run this if you have TreeSummarizedExperiment installed)
if (FALSE) { # \dontrun{
library("TreeSummarizedExperiment")
example("TreeSummarizedExperiment")
assayNames(tse) <- "counts"
emuRes <- emuFit(Y = tse, formula = ~ condition, assay_name = "counts", 
                 test_kj = data.frame(k = 2, j = 1), tolerance = 0.01)
 # here we set large tolerances for the example to run quickly, 
 # but we recommend smaller tolerances in practice
} # }
```
