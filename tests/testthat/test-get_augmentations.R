test_that("in saturated case, augmentations reduce to haldane correction", {
    Y <- matrix(1:4,nrow = 2)
    X <- cbind(1,c(0,1))
    J <- ncol(Y)
    n <- nrow(Y)


    ml_fit <- emuFit_micro(X,Y ,
                           constraint_fn = rep(list(function(x) x[2]), 2), 
                           tolerance = 1e-5,
                           verbose = FALSE)

    X_cup <- X_cup_from_X_fast(X,J)
    G <- get_G_for_augmentations_fast(X = X,J = J,n = n,
                                 X_cup = X_cup)
    expect_true(max(abs(get_augmentations(X = X,
                                          G = G,
                                          Y = Y,
                                          B = ml_fit) - 0.5)) < 1e-10)

    #correction should not depend on B in this case -- check this too
    expect_true(max(abs(get_augmentations(X = X,
                                          G = G,
                                          Y = Y,
                                          B = 0*ml_fit) - 0.5)) < 1e-10)

  })


# # test timing (remove this later)
# 
# meta <- readRDS("../crystals_analysis/abx_recency_taxa_count_metadata_dummies.rds")
# abundance <- readRDS("../crystals_analysis/abx_recency_taxa_count_dummies.rds")
# abundance <- abundance[, 1:300]
# design <- make_design_matrix(formula = ~ abx_recency_three + abx_recency_twelve + 
#                                center_mean_age + sex, data = meta)
# X_cup <- X_cup_from_X_fast(design, ncol(abundance))
# G <- get_G_for_augmentations_fast(design, ncol(abundance), nrow(abundance), X_cup)
# debugonce(get_augmentations)
# j_ref <- get_j_ref(abundance)
# Y_augmented <- as.matrix(abundance) + 1e-3*mean(as.numeric(unlist(abundance)))
# fitted_model <- emuFit_micro(design,
#                              Y_augmented,
#                              j_ref = j_ref, maxit = 5)
# start1 <- proc.time()
# aug1 <- get_augmentations(design, G, as.matrix(abundance), fitted_model)
# end1 <- proc.time() - start1
# 
# start2 <- proc.time() 
# aug2 <- get_augmentations_fast(design, G, as.matrix(abundance), fitted_model)
# end2 <- proc.time() - start2
# 
# all.equal(aug1, aug2)
# 
# 
# 
# 
# meta <- readRDS("../crystals_analysis/abx_recency_taxa_count_metadata_dummies.rds")
# meta <- meta[1:200, ]
# abundance <- readRDS("../crystals_analysis/abx_recency_taxa_count_dummies.rds")
# abundance <- as.matrix(abundance[1:200, ])
# design <- make_design_matrix(formula = ~ abx_recency_three + abx_recency_twelve + 
#                                center_mean_age + sex + center_mean_bmi + 
#                                center_mean_gfr + meds_blood_sugar + 
#                                meds_cholesterol, data = meta)
# X_cup <- X_cup_from_X_fast(design, ncol(abundance))
# G <- get_G_for_augmentations_fast(design, ncol(abundance), nrow(abundance), X_cup)
# j_ref <- get_j_ref(abundance)
# Y_augmented <- as.matrix(abundance) + 1e-3*mean(as.numeric(unlist(abundance)))
# new_fitted_model <- emuFit_micro(design,
#                              Y_augmented,
#                              j_ref = j_ref, maxit = 5)
# start1 <- proc.time()
# aug1 <- get_augmentations(design, G, as.matrix(abundance), new_fitted_model)
# end1 <- proc.time() - start1
# 
# start2 <- proc.time() 
# aug2 <- get_augmentations_fast(design, G, as.matrix(abundance), new_fitted_model)
# end2 <- proc.time() - start2
# 
# all.equal(aug1, aug2)
# end1
# end2
