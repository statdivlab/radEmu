# Calculate the derivative to the pseudo-Huber smoothed median

Calculate the derivative to the pseudo-Huber smoothed median

## Usage

``` r
dpseudohuber_median_dx(x, d = 0.1, na.rm = FALSE)

dpsuedohuber_median_dx(x, d = 0.1, na.rm = FALSE)
```

## Arguments

- x:

  A vector to calculate the derivative of the pseudo-Huber smoothed
  median for.

- d:

  Smoothing parameter, by default set to `0.1`. As `d` approaches `0`
  the pseudo-Huber median function approaches the median and as `d`
  approaches infinity this function approaches the mean.

- na.rm:

  Passed to `pseudohuber_median`, default is FALSE, if FALSE then when
  `x` includes at least one NA value then NA is returned, if TRUE then
  when `x` includes at least one NA value then that value is removed and
  the pseudo-Huber median is computed without it.

## Value

The derivative of the calculated pseudo-Huber smoothed median over `x`
with smoothing parameter `d`.

## Examples

``` r
dpseudohuber_median_dx(x = rnorm(10), d = 0.1)
#>  [1] 3.247388e-03 3.512036e-05 2.745777e-02 6.179727e-04 1.464490e-04
#>  [6] 9.309298e-05 3.532373e-01 3.601095e-01 2.456994e-01 9.355974e-03
```
