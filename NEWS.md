# radEmu 2.1.1.0

This is a minor release that exports the functions `pseudohuber_median()` and `dpseudohuber_median_dx()`, which calculate the pseudo-Huber smoothed median and its derivative. These are the default constraint functions used in `radEmu`. 

# radEmu 2.0.0.0

This is a major release that speeds up score tests, and forces the user to clarify that they wish to perform score tests. It makes the default behaviour faster, but is not backwards compatible.

## Breaking changes

* The argument `test_kj` is now required for `emuFit()` when `run_score_tests = TRUE` (the default). Previous default behavior was to run score tests for every parameter, which can be very time consuming(and can easily be parallelized). This change forces the user to explicitly state what tests they would like to run, significantly decreasing unnecessary computation. 

## Additional changes

* We have also streamlined estimation under the null, leading to improved convergence and reduced computation. 