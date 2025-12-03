# Data simulation function

Function to simulate data for simulations in Clausen & Willis (2024) and
for the cluster vignette

## Usage

``` r
simulate_data(
  n,
  J,
  b0 = NULL,
  b1 = NULL,
  distn,
  zinb_size = NULL,
  zinb_zero_prop = NULL,
  mean_z,
  X = NULL,
  B = NULL,
  cluster = NULL
)
```

## Arguments

- n:

  Number of samples

- J:

  Number of categories

- b0:

  Intercept parameter vector

- b1:

  Covariate paramter vector

- distn:

  Distribution to simulate from, either "Poisson" or "ZINB"

- zinb_size:

  Size parameter for negative binomial draw for ZINB data

- zinb_zero_prop:

  Proportion of zeros for ZINB data

- mean_z:

  Parameter controlling the mean of the sample-specific effects.

- X:

  Optional design matrix, this must have two columns and n rows.

- B:

  Optional B matrix, if p is not equal to 2

- cluster:

  Optional cluster vector, this must have n elements.

## Value

`Y`. A `n times J` dimension matrix of simulated response counts.
