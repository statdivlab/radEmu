#' A fitted GLM model for the power of radEmu to reject a truly false null hypothesis
#'
#' 100 simulations were drawn with n = 10, 30, 50 and 100 total samples; J = 250 and 500 taxa; 
#' a single categorical covariate (p = 2); effect sizes beta11 (the parameter of interest) from 
#' 0, 0.5, 1, ... 2.5. beta0's range from -3 to 3, and the other beta1's range from -1 to 1
#' with no correlation between the beta0s (which, roughly speaking, control the relative abundance of the taxa) and
#' the beta1s (which control the difference in abundance between the two covariate groups). 
#' Counts were drawn from a zero-inflated negative binomial model with size parameter 5, zero-inflation 
#' probability of 0.5 and average z's around log(50). 
#' A model for the probability of rejecting the null hypothesis of beta11 = 0 was fit. Model fitting was guided by 
#' plotting the log odds of rejection, where effect modification between n and beta11 was observed. 
#' This model may be useful for power calculations in future, though as with any simulation, 
#' its generalizability is limited to similar data generating processes. Simulation code that can be generalized is available at
#' https://github.com/statdivlab/radEmu_supplementary under fig-power/power_simulations.R
#'
#' @format A GLM object.
#' \describe{
#' \item{power_model}{A GLM object modelling the odds of rejecting the null hypothesis at a given sample size, number of taxa, and effect size}
#' }
#' 
#' @references Wirbel, J et al. (2023). \emph{Meta-analysis of fecal metagenomes reveals global microbial signatures that are specific for colorectal cancer}. Nature Medicine, 25, 679–689. <doi: 10.1038/s41591-019-0406-6>.
"power_model"