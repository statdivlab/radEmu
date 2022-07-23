
do_one_sim <- function(n_per_group,
                       J,
                       delta,
                       b1_magnitude,
                       cor_structs,
                       seed,
                       nboot = 250,
                       verbose = FALSE

){

  b1_pre <- c(1,rep(0,J/4 - 1),
              -1,rep(0,J/4 - 1),
              rep(0,J/4 - 1), 1,
              rep(0,J/4 - 1),-1)

  set.seed(seed)

  Y <- generate_sim_data(n_per_group = n_per_group,
                    b0 = c(rep(-9,J/10),seq(-5,-2,length.out = 0.8*J),rep(2,J/10)),
                    b1 = b1_pre*b1_magnitude,
                    delta = delta,
                    z_mean = 4,
                    nb_size = .125,
                    cor_structs = cor_structs)

  absent <- apply(Y,2,function(x) sum(x) ==0)
  which_absent <- which(absent)
  X <- cbind(1,rep(c(0,1),each = n_per_group))

  message("Removing absent taxa")
  Y <- Y[,!absent]

  ml_fit <- emuFit(X = X, Y = Y,method = "ML",
                   verbose = verbose)
  ml_reweighted_fit <- emuFit(X = X,
                              Y = Y, method = "ML",
                              B = ml_fit$B,
                              reweight = TRUE,
                              prefit = FALSE,
                              verbose = verbose)

  fl_fit <- emuFit(X = X, Y = Y,
                   method = "FL",
                   verbose = verbose)

  fl_reweighted_fit <- emuFit(X = X,
                              Y = Y,
                              B = fl_fit$B,
                              method = "FL",
                              reweight = TRUE,
                              prefit = FALSE,
                              verbose = verbose)

  ml_cis <- emuCI(ml_fit,nboot = nboot,
                  rows_of_interest = 2,
                  verbose = verbose)

  ml_reweighted_cis <- emuCI(ml_reweighted_fit,nboot = nboot,
                  rows_of_interest = 2,
                  verbose = verbose)

  fl_cis <- emuCI(fl_fit,nboot = nboot,
                  rows_of_interest = 2,
                  verbose = verbose)

  fl_reweighted_cis <- emuCI(fl_reweighted_fit,nboot = nboot,
                  rows_of_interest = 2,
                  verbose = verbose)

  n_missing <- length(which_absent)

  message("Reinserting rows for absent taxa")

  ml_cis <- sim_insert_absent(ci_object = ml_cis,
                    which_absent = which_absent,
                    J = J)

  ml_reweighted_cis <- sim_insert_absent(ci_object = ml_reweighted_cis,
                                         which_absent = which_absent,
                                         J = J)

  fl_cis <- sim_insert_absent(ci_object = fl_cis,
                              which_absent = which_absent,
                              J = J)

  fl_reweighted_cis <- sim_insert_absent(ci_object = fl_reweighted_cis,
                                          which_absent = which_absent,
                                         J = J)

  message("Collating results")
  ml_cis$estimator <- "ml"
  ml_cis$weighted <- "no"
  ml_reweighted_cis$estimator <- "ml"
  ml_reweighted_cis$weighted <- "yes"
  fl_cis$estimator <- "fl"
  fl_cis$weighted <- "no"
  fl_reweighted_cis$estimator <- "fl"
  fl_reweighted_cis$weighted <- "yes"

  return(rbind(ml_cis,
               ml_reweighted_cis,
               fl_cis,
               fl_reweighted_cis))

}
