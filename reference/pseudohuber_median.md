# Calculate the pseudo-Huber smoothed median

Calculate the pseudo-Huber smoothed median using a quadratic
approximation to the pseudo-Huber criterion detailed in supplement to
Clausen et al. (2024).

## Usage

``` r
pseudohuber_median(x, d = 0.1, tolerance = 1e-08, na.rm = FALSE)

psuedohuber_median(x, d = 0.1, tolerance = 1e-08, na.rm = FALSE)
```

## Arguments

- x:

  A vector to calculate the pseudo-Huber smoothed median for.

- d:

  Smoothing parameter, by default set to `0.1`. As `d` approaches `0`
  this function approaches the median and as `d` approaches infinity
  this function approaches the mean.

- tolerance:

  Tolerance used to determine convergence in the algorithm used to
  calculate this function value.

- na.rm:

  Default is FALSE, if FALSE then when `x` includes at least one NA
  value then NA is returned, if TRUE then when `x` includes at least one
  NA value then that value is removed and the pseudo-Huber median is
  computed without it.

## Value

The calculated pseudo-Huber smoothed median over `x` with smoothing
parameter `d`.

## Examples

``` r
pseudohuber_median(x = rnorm(10), d = 0.1)
#> [1] -0.07817928
```
