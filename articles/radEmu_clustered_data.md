# Using radEmu on clustered data

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

## Introduction

In this vignette we will explore the use of `radEmu` with clustered
data. Data with cluster dependence is a common phenomenon in microbiome
studies. This type of dependence can come from experimental factors such
as shared cages or tanks for study animals. It can also arise from
repeated measurements in longitudinal studies.

Sometimes people deal with clustered data via random effects. You might
have seen the syntax `+ (1 | cluster)`. We (the `radEmu` developers)
don’t want to make strong assumptions (such as the random effects being
normally distributed), so we handle cluster dependence using a GEE
framework. We like this approach because it’s robust to many forms of
model misspecification, and for this reason, we find it superior to
random effects.

When dependence is not accounted for, statistical inference in `radEmu`
(and most statistical methods!) assume that all samples are independent.
If you have cluster dependence but assume independence, you’ll have
anti-conservative inference (i.e., p-values that are smaller than they
should be). **Therefore, we strongly recommend adjusting for cluster
dependence if it arose in your sample collection!**

Note that cluster dependence won’t change your estimates
(${\widehat{\beta}}_{j}$’s), but it will (most likely) change your
p-values.

Luckily, we have tools to account for cluster dependence implemented in
`radEmu`! This argument is only implemented from `radEmu` v1.2.0.0
forward, so if you are having trouble using the `cluster` argument,
check that you reinstalled `radEmu` recently.

TLDR; the basic syntax is as follows. `my_clusters` is a vector of
length $n$ (your number of samples), with observations from the same
cluster having the same value in `my_clusters`
(e.g. `my_clusters = c(1, 1, 2, 2, 3, 4)`).

    emuFit(formula = ~ covariate,
           data = my_data_frame, 
           Y = my_microbial_abundances,
           cluster = my_clusters)

Fun fact! We implemented this functionality because of user requests.
Therefore, if there’s something that you’d like to see that you don’t
see, [let us know](https://github.com/statdivlab/radEmu/issues) and
we’ll see what we can do!

## Generating data with cluster dependence

To start, let’s generate a toy example (10 categories, 60 samples) in
which there is cluster dependence within our data. The way we simulate
data isn’t important; it’s just an illustration. **radEmu can handle
more taxa and samples,** we just did this so that the vignette builds
quickly.

``` r
J <- 10; n <- 60
# generate design matrix 
set.seed(10)
X <- cbind(1, rnorm(n))
cov_dat <- data.frame(cov = X[, 2])
# cluster membership 
cluster <- rep(1:4, each = n/4)
cluster_named <- paste("cage", cluster, sep = "")
cov_dat$cluster <- cluster
# intercepts for each category
b0 <- rnorm(J)
# coefficients for X1 for each category 
b1 <- seq(1, 5, length.out = J)
# mean center the coefficients
b1 <- b1 - mean(b1)
# set the coefficient for the 3rd category to 4 (why not!?)
# Note that because of the constraint, we're only able to estimate
# b1 - mean(b1), which is ~3.9.  
b1[3] <- 4
# generate B coefficient matrix 
b <- rbind(b0, b1)

# simulate data according to a zero-inflated negative binomial distribution
# the mean model used to simulate this data takes into account the cluster membership
Y <- radEmu:::simulate_data(n = n,
                            J = J,
                            b0 = b0,
                            b1 = b1,
                            distn = "ZINB",
                            zinb_size = 10,
                            zinb_zero_prop = 0.3,
                            mean_z = 5,
                            X = X,
                            cluster = cluster)
```

Let’s just pause to look at the elements of `cluster`:

``` r
table(cluster_named)
#> cluster_named
#> cage1 cage2 cage3 cage4 
#>    15    15    15    15
```

So all of the observations from the first cage (or person, tank,
whatever…) have `cage1` in their corresponding `cluster_named` variable.

Let’s fit a model to this data! We know that the log-fold difference
between in the abundance of category 3 when comparing samples that
differ by 1 unit in $X$ is 3.5, so if we have good power, we will reject
the null that $\beta_{X_{1},3} = 0$. We fit a model including cluster
dependence as follows:

``` r
ef_cluster <- emuFit(formula = ~ cov,
                     data = cov_dat, 
                     Y = Y,
                     cluster = cluster_named, 
                     test_kj = data.frame(k = 2, j = 3))
#> Row names are missing from the response matrix Y. We will assume the rows are in the same order as in the covariate matrix X. You are responsible for ensuring the order of your observations is the same in both matrices.
```

You can check out the full object, but our estimate is 4.2 and a p-value
for testing that this parameter equals zero is 0.03. Not too shabby,
especially considering that about half of our observations are zero, and
we have a lot of noise in our data (arising from a negative binomial
simulation scheme).

Let’s also compare that to a situation where we mistakenly ignore
clustering. In this case, we expect to have a smaller p-value, because
we are saying that we have more independent observations.

``` r
ef_no_cluster <- emuFit(formula = ~ cov,
                        data = cov_dat, 
                        Y = Y,
                        test_kj = data.frame(k = 2, j = 3))
#> Row names are missing from the response matrix Y. We will assume the rows are in the same order as in the covariate matrix X. You are responsible for ensuring the order of your observations is the same in both matrices.
```

When we ignore clustering, we get an estimate of the log fold difference
in category 3 across values of our covariate of 4.2 and a p-value of
0.005. So we can see that our estimates are the same whether or not we
account for cluster, but our p-values are different (because we are
pretending that we have more evidence in the absence of clustering).
