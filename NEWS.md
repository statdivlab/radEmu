# radEmu 2.0.0.0

This is a major release that updates the argument `test_kj` in the function `emuFit()` in a way that is not backwards compatible.

## Breaking changes

* When `emuFit()` is run with the default of `run_score_tests = TRUE`, the argument `test_kj` is now required instead of optional. In previous versions when`run_score_tests = TRUE` and `test_kj` was not included, score tests were run for all parameters in the model except for intercept parameters. As these score tests can be computationally expensive, we want to make sure that the user only runs tests for parameters that they want to include in the analysis. 
