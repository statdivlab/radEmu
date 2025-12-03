# Parallelizing computation for score tests with radEmu

First, we will install `radEmu`, if we haven’t already.

``` r
# if (!require("remotes", quietly = TRUE))
#     install.packages("remotes")
# 
# remotes::install_github("statdivlab/radEmu")
```

Next, we can load `radEmu` as well as the `tidyverse` package suite.

``` r
library(magrittr)
library(dplyr)
library(ggplot2)
library(stringr)
library(radEmu)
```

Finally, we will load the package `parallel`.

``` r
library(parallel)
```

## Introduction

In this vignette we will introduce parallel computing in order to do
more efficient computation for score tests. Our recommendation for
testing statistical hypotheses with small to moderate sample sizes with
`radEmu` is to run a robust score test. While this test performs well
(we like that it controls Type I error rate even in small samples!), it
takes some time to run, because we need to fit the model under each null
hypothesis. For differential abundance analysis, we often want to run a
hypothesis test for each category (taxon, gene, etc) that we care about,
so this adds up quickly. In order to improve computational efficiency,
we can run these score tests in parallel using the `parallel` R package.
This will let us take advantage of having additional cores on personal
computers or computing clusters. Note that we will be using the
`mclapply` function from the `parallel` package, which works on Mac and
Linux machines, but does not on Windows. If you are using a Windows
machine and would like a vignette about parallel computing on Windows,
please let us know by opening an issue.

We recommend that before working through this vignette, you start the an
introduction to `radEmu` in “intro_radEmu.Rmd.”

## Setting up our radEmu model and running a single score test

We’ll use the same Wirbel et al. data as in the introduction vignette.
Recall that the [dataset published by Wirbel et
al. (2019)](https://www.nature.com/articles/s41591-019-0406-6) is from a
meta-analysis of case-controls comparing participants with and without
colorectal cancer.

``` r
# load in sample data
data("wirbel_sample")
# set group to be a factor with levels CTR for control and CRC for cancer
wirbel_sample$Group <- factor(wirbel_sample$Group, levels = c("CTR","CRC"))
# load in abundance data 
data("wirbel_otu")
# save mOTU names
mOTU_names <- colnames(wirbel_otu)
# consider taxa in the following genera
chosen_genera <- c("Eubacterium", "Faecalibacterium", "Fusobacterium", "Porphyromonas")
# get taxonomy information from mOTU names
mOTU_name_df <- data.frame(name = mOTU_names) %>% 
  mutate(base_name = stringr::str_remove(mOTU_names, "unknown ") %>%
                      stringr::str_remove("uncultured ")) %>%
  mutate(genus_name = stringr::word(base_name, 1))
# restrict to names in chosen genera
restricted_mOTU_names <- mOTU_name_df %>%
  filter(genus_name %in% chosen_genera) %>%
  pull(name)
# pull out observations from a chinese study within the meta-analysis
ch_study_obs <- which(wirbel_sample$Country %in% c("CHI"))
# make count matrix for chosen samples and genera
small_Y <- wirbel_otu[ch_study_obs, restricted_mOTU_names]
# check for samples with only zero counts
sum(rowSums(small_Y) == 0) # no samples have a count sum of 0 
#> [1] 0
# check for genera with only zero counts
sum(colSums(small_Y) == 0) # one category has a count sum of 0
#> [1] 1
# remove the one genus with only zero counts
category_to_rm <- which(colSums(small_Y) == 0)
small_Y <- small_Y[, -category_to_rm]
```

Now that we’ve processed our data, we can fit the `radEmu` model. Here
we just want to get estimates for our parameters and their standard
errors, but we will avoid running score tests by setting
`run_score_tests = FALSE`.

``` r
ch_fit <- emuFit(formula = ~ Group, 
                 data = wirbel_sample[ch_study_obs, ],
                 Y = small_Y,
                 run_score_tests = FALSE) 
```

In the introduction vignette we found that a meta-mOTU “unknown
Eubacterium \[meta_mOTU_v2_7116\]” assigned to Eubacteria has a much
higher ratio of abundance (comparing CRC group to control) than is
typical across the mOTUs we included in this analysis, based on the
parameter estimates in `ch_fit`. We can run a robust score test to test
whether the differential abundance of this mOTU between cases and
controls is significantly different from the differential abundance of a
typical mOTU in our analysis.

In order to run a single robust score test, we will set
`run_score_tests = TRUE` and include the argument `test_kj`. Instead of
re-estimating parameters in our model, we will provide `ch_fit` to the
argument `fitted_model` and set `refit = FALSE`.

``` r
mOTU_to_test <- which(str_detect(restricted_mOTU_names, "7116"))
ch_fit$B %>% rownames
#> [1] "(Intercept)" "GroupCRC"
covariate_to_test <- which("GroupCRC" == ch_fit$B %>% rownames)
robust_score <- emuFit(formula = ~ Group,
                       data = wirbel_sample[ch_study_obs, ],
                       fitted_model = ch_fit,
                       refit = FALSE,
                       test_kj = data.frame(k = covariate_to_test, 
                                            j = mOTU_to_test), 
                       Y = small_Y)
robust_score$coef$pval[mOTU_to_test]
#> [1] 0.3017219
```

Now, we can see that it took a little while to run our robust score
test. If we investigate the coefficient table in our `robust_score`
output, we can see that we have a p-value of
`R robust_score$coef$pval[mOTU_to_test]` from our test.

## Running robust score tests in parallel

Now, let’s run some tests in parallel. We will be parallelizing our code
over $j$ in the argument `test_kj`. We will assume that you have one
covariate that you want to test, corresponding with a specific column of
your design matrix $k$. However, if you want to run tests for multiple
columns of your design matrix, then you can parallelize over pairs of
$k$ and $j$ in the argument $test_{k}j$.

Let’s say that I want to run score tests for the first five mOTUs in our
dataset. First, I need to check the cores on my computer to see a
reasonable amount of cores to parallelize over. I tend to use one fewer
core than the number of cores that I have available.

``` r
ncores <- parallel::detectCores() - 1
ncores
#> [1] 3

# if running this vignette in automatic R-CMD-check, reduce cores to 2
if (identical(Sys.getenv("_R_CHECK_LIMIT_CORES_"), "TRUE")) {
  ncores <- min(2, ncores)
}
```

Next, I will write a function that will be called in parallel. This
function will fit the model under the null and calculate the robust
score test statistics. Note that the output of this function is an
`emuFit` object.

``` r
emuTest <- function(category) {
  score_res <- emuFit(formula = ~ Group,
                       data = wirbel_sample[ch_study_obs, ],
                       fitted_model = ch_fit,
                       refit = FALSE,
                       test_kj = data.frame(k = covariate_to_test, 
                                            j = category), 
                       Y = small_Y)
  return(score_res)
}
```

Now, we can run our score tests in parallel. We’ll just do the first
five. It may take a minute or so, depending on your machine.

``` r
if (.Platform$OS.type != "windows" & !identical(Sys.getenv("GITHUB_ACTIONS"), "true")) {
  # run if we are on a Mac or Linux machine
  score_res <- mclapply(1:5,
                      emuTest,
                      mc.cores = ncores)
} else {
  # don't run if we are on a Windows machine, or if testing with GitHub actions
  score_res <- NULL
}
```

Now, we can see that this barely took more time than running a single
score test, because we were able to parallelize over cores (on my
laptop, I’m using more than five cores, so I can run all five tests at
the same time). Each p-value can be pulled out of the list as follows:

``` r
if (!is.null(score_res)) {
  c(score_res[[1]]$coef$pval[1], ## robust score test p-value for the first taxon
  score_res[[2]]$coef$pval[2]) ## robust score test p-value for the second taxon
}
```

To help organise this information, we can make a coefficient matrix that
combines the information from each component in our list:

``` r
if (!is.null(score_res)) {
  full_score <- sapply(1:length(score_res), 
                       function(x) score_res[[x]]$coef$score_stat[x])
  full_pval <- sapply(1:length(score_res), 
                      function(x) score_res[[x]]$coef$pval[x])
  full_coef <- ch_fit$coef %>%
    dplyr::select(-score_stat, -pval) %>%
    filter(category_num %in% 1:5) %>%
    mutate(score_stat = full_score,
           pval = full_pval)
  full_coef
}
```

The column containing our p-values is called `pval`.

Happy testing!
