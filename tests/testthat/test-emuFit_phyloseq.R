library(radEmu)

test_that("emuFit works with a phyloseq object", {
  if (requireNamespace("phyloseq", quietly = TRUE)) {
    # make phyloseq objects from data
    data(wirbel_sample)
    data(wirbel_otu)
    data(wirbel_taxonomy)
    wirbel_phylo <- phyloseq::phyloseq(phyloseq::sample_data(wirbel_sample),
                                       phyloseq::otu_table(wirbel_otu, taxa_are_rows = FALSE),
                                       phyloseq::tax_table(wirbel_taxonomy))
    wirbel_small <- phyloseq::tax_glom(wirbel_phylo, "phylum") 
    wirbel_smaller <- phyloseq::prune_samples(wirbel_sample$Country == "FRA" &
                                                wirbel_sample$Gender == "F", 
                                              wirbel_small)
    fit <- emuFit(wirbel_small, formula = ~ Group, run_score_tests = FALSE, tolerance = 0.01)
    expect_true(is.matrix(fit$B))
    
    wirbel_transpose <- wirbel_small
    phyloseq::otu_table(wirbel_transpose) <- t(phyloseq::otu_table(wirbel_transpose))
    fit_transpose <- emuFit(wirbel_transpose, formula = ~ Group, run_score_tests = FALSE, tolerance = 0.01)
    expect_true(all.equal(fit$coef, fit_transpose$coef))
    
  } else {
    expect_error(stop("You don't have phyloseq installed"))
  }
})
