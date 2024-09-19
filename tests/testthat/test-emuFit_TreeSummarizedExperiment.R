library(radEmu)

test_that("emuFit works with a TreeSummarizedExperiment object", {
  if (requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
    # make TreeSummarizedExperiment object from data
    data(wirbel_sample)
    data(wirbel_otu)
    data(wirbel_taxonomy)
    
    sub_samples <- wirbel_sample$Country == "FRA" & wirbel_sample$Gender == "F"
    sub_taxa <- colSums(wirbel_otu[sub_samples, , drop = FALSE]) > 0
    
    wirbel_sample_sub <- wirbel_sample[sub_samples, ]
    wirbel_otu_sub <- wirbel_otu[sub_samples, sub_taxa]
    wirbel_taxonomy_sub <- wirbel_taxonomy[sub_taxa, ]
    
    Y <- TreeSummarizedExperiment::TreeSummarizedExperiment(
      assays = list("counts" = t(wirbel_otu_sub)),
      rowData = wirbel_taxonomy_sub,
      colData = wirbel_sample_sub)
    
    # fit model using TreeSummarizedExperiment
    fit <- emuFit(Y = Y,
                  formula = ~ Group,
                  assay_name = "counts",
                  run_score_tests = FALSE, tolerance = 0.01)
    
    # fit model using data.frames directly
    fit2 <- emuFit(Y = wirbel_otu_sub,
                   X = model.matrix(object = ~ Group, data = wirbel_sample_sub),
                   run_score_tests = FALSE, tolerance = 0.01)
    
    # confirm the results match when data are extracted from TreeSummarizedExperiment
    # or when data.frames are used directly
    expect_true(all.equal(fit$coef, fit2$coef))
    
  } else {
    expect_error(stop("You don't have TreeSummarizedExperiment installed."))
  }
})
