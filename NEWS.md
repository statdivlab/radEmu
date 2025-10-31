# radEmu 2.2.0.0

This is a minor release that changes the algorithm used to estimate parameters under the null hypothesis in score tests depending on the identifiability constraint that is used. All functions are backwards compatible, but there are major changes in the back-end, which may lead to small changes in results of analyses compared to those that were run before this update. 

## Changes

A different algorithm is used by default to estimate parameters under the null hypothesis. This algorithm is typically faster and more likely to converge to a solution.

Additionally, a `control` argument list is added to `emuFit()`, which includes many arguments that relate to algorithms and convergence. This is not a breaking change because using these arguments outside of the `control` list will not lead to an error. The goal of this is to clarify documentation. 

# radEmu 2.1.1.0

This is a minor release that exports the functions `pseudohuber_median()` and `dpseudohuber_median_dx()`, which calculate the pseudo-Huber smoothed median and its derivative. These are the default constraint functions used in `radEmu`. 

# radEmu 2.0.0.0

This is a major release that speeds up score tests, and forces the user to clarify that they wish to perform score tests. It makes the default behaviour faster, but is not backwards compatible.

## Breaking changes

* The argument `test_kj` is now required for `emuFit()` when `run_score_tests = TRUE` (the default). Previous default behavior was to run score tests for every parameter, which can be very time consuming(and can easily be parallelized). This change forces the user to explicitly state what tests they would like to run, significantly decreasing unnecessary computation. 

## Additional changes

* We have also streamlined estimation under the null, leading to improved convergence and reduced computation. 