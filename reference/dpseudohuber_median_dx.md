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
#>  [1] 2.328500e-04 3.212551e-03 3.513189e-05 2.694470e-02 6.136539e-04
#>  [6] 1.457170e-04 9.319663e-05 3.560505e-01 3.625244e-01 2.501474e-01
```
