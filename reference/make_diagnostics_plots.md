# makes plots to investigate convergence of estimation under null hypothesis

makes plots to investigate convergence of estimation under null
hypothesis

## Usage

``` r
make_diagnostics_plots(diagnostic_df)
```

## Arguments

- diagnostic_df:

  Dataframe with relevant information

## Value

Plots with diagnostic information. If the fisher scoring algorithm has
been fit, this will be two plots showing how the log likelihood and test
statistic change over iterations. If the augmented lagrangian algorithm
has been fit, this will also include how the constraint gap and maximum
changing element of B change over iterations.
