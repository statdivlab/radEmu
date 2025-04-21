# radEmu 3.0.0.0

This is a major release that changes the algorithm used to estimate parameters under the null hypothesis in score tests depending on the identifiability constraint that is used. All functions are backwards compatible, but there are major changes in the backend, which may lead to small changes in results of analyses compared to those that were run before this update. 

# radEmu 2.0.0.0

This is a major release that speeds up score tests, and forces the user to clarify that they wish to perform score tests. It makes the default behaviour faster, but is not backwards compatible.

## Breaking changes

* The argument `test_kj` is now required for `emuFit()` when `run_score_tests = TRUE` (the default). Previous default behavior was to run score tests for every parameter, which can be very time consuming(and can easily be parallelized). This change forces the user to explicitly state what tests they would like to run, significantly decreasing unnecessary computation. 

## Additional changes

* We have also streamlined estimation under the null, leading to improved convergence and reduced computation. 