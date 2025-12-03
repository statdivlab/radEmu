# Generates the design matrix for a given `formula` and `data` object

Generates the design matrix for a given `formula` and `data` object

## Usage

``` r
make_design_matrix(Y = NULL, formula, data = NULL, assay_name = NULL)
```

## Arguments

- Y:

  Optionally, a `phyloseq` or `TreeSummarizedExperiment` object
  containing an otu table and sample data.

- formula:

  a one-sided formula specifying the form of the mean model to be fit

- data:

  an n x p data frame containing variables given in `formula` (required
  unless `Y` is included as a `phyloseq` or `TreeSummarizedExperiment`
  object)

- assay_name:

  Optionally, a string containing the desired assay name within a
  `TreeSummarizedExperiment` object. This is only required if `data` is
  a `TreeSummarizedExperiment` object, otherwise this argument does
  nothing and can be ignored.

## Value

returns design matrix `X`.
