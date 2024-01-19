#' Covariates from Wirbel et al. meta-analysis of fecal metagenomes.
#'
#' A data frame with covariates. 
#'
#' @format A data frame object with sample data, 566 observations of 14 covariates.
#' \describe{
#' \item{wirbel_sample}{sample data with the following covariates:
#' \itemize{
#' \item \code{Sample_ID}, id of specific sample
#' \item \code{External_ID}, id of specific sample from original study
#' \item \code{Age}, age of subject from which sample was taken
#' \item \code{Gender}, gender of subject from which sample was taken
#' \item \code{BMI}, BMI of subject from which sample was taken 
#' \item \code{Country}, country of study
#' \item \code{Study}, study that sample is from 
#' \item \code{Group}, CRC for colorectal cancer or CTR for control 
#' \item \code{Library_Size}, library size of sample
#' \item \code{Age_spline.1}, value of first coordinate of age spline
#' \item \code{Age_spline.2}, value of second coordinate of age spline
#' \item \code{BMI_spline.1}, value of first coordinate of BMI spline
#' \item \code{BMI_spline.2}, value of second coordinate of age spline
#' \item \code{Sampling}, whether sampling happened before or after colonoscopy
#' }}
#' }
#' @references Wirbel, J et al. (2019). \emph{Meta-analysis of fecal metagenomes reveals global microbial signatures that are specific for colorectal cancer}. Nature Medicine, 25, 679–689. <doi: 10.1038/s41591-019-0406-6>.
"wirbel_sample"