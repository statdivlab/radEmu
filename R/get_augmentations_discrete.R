#get data augmentations for Firth penalized estimation in closed form, for
#designs with only discrete covariates (X has exactly p distinct rows)
#equal to get_augmentations(X, G, Y, B) for such designs: with z restricted as
#in get_augmentations, the hat value for observation (i, j) is
#p_ij + (1 - p_ij) * ||Y_i||_1 / N_g, where p_ij is the fitted proportion of
#category j for sample i and N_g is the total of Y over samples sharing
#sample i's covariate pattern g
get_augmentations_discrete <- function(X,
                                       Y,
                                       B,
                                       groups = discrete_groups(X)){

  #fitted proportions, one row per covariate pattern
  log_props <- groups$distinct_X %*% B
  props <- exp(log_props - apply(log_props, 1, max))
  props <- props / rowSums(props)

  Y_rowsums <- rowSums(Y)
  group_totals <- as.numeric(rowsum(Y_rowsums, groups$group, reorder = TRUE))
  w <- Y_rowsums / group_totals[groups$group]

  #hat values are props * (1 - w) + w; augmentations are half the hat values
  (props[groups$group, , drop = FALSE] * (1 - w) + w) / 2
}
