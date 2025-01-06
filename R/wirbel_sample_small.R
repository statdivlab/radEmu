#' Covariates from Wirbel et al. meta-analysis of fecal metagenomes, for a subset of samples.
#'
#' A data frame with covariates. 
#'
#' @format A data frame object with sample data, 566 observations of 14 covariates.
#' \describe{
#' \item{Sample_ID}{id of specific sample}
#' \item{External_ID}{id of specific sample from original study}
#' \item{Age}{age of subject from which sample was taken}
#' \item{Gender}{gender of subject from which sample was taken}
#' \item{BMI}{BMI of subject from which sample was taken}
#' \item{Country}{country of study}
#' \item{Study}{study that sample is from} 
#' \item{Group}{CRC for colorectal cancer or CTR for control} 
#' \item{Library_Size}{library size of sample}
#' \item{Age_spline.1}{value of first coordinate of age spline}
#' \item{Age_spline.2}{value of second coordinate of age spline}
#' \item{BMI_spline.1}{value of first coordinate of BMI spline}
#' \item{BMI_spline.2}{value of second coordinate of age spline}
#' \item{Sampling}{whether sampling happened before or after colonoscopy}
#' }
#' @references Wirbel, J et al. (2019). \emph{Meta-analysis of fecal metagenomes reveals global microbial signatures that are specific for colorectal cancer}. Nature Medicine, 25, 679–689. <doi: 10.1038/s41591-019-0406-6>.
"wirbel_sample_small"