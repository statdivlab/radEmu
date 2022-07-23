test_that("Empirical means in simulated model are consistent with theoretical values", {

  J <- 50
  n_per_group <- 50
  set.seed(0)
  cor_structs <-
    generate_sim_setting(J,
                       n_correlations = 20,
                       max_no_correlated = 10)

  delta <- rnorm(J,sd = 3)


 draws <- lapply(1:100,
        function(k)
          generate_sim_data(n_per_group = 50,
                    b0 = rep(-1,J),
                    b1 = rep(0,J),
                    delta = delta,
                    z_mean = 6,
                    nb_size = 1,
                    cor_structs = cor_structs,
                    constant_z = TRUE))

  means <- matrix(nrow = 2*n_per_group,
                  ncol = J)
 for(i in 1:(2*n_per_group)){
   for(j in 1:J){
     means[i,j] <- mean(sapply(1:100,
                               function(k) draws[[k]][i,j]))
   }
 }

  expect_true(abs(lm(log(apply(means,2,mean))~delta)$coef[2] - 1)< 0.05)
})

#
test_that("Under null simulation settings, observed counts are reasonable", {

  J <- 50
  n_per_group <- 25
  set.seed(0)
  cor_structs <-
    generate_sim_setting(J,
                         n_correlations = 20,
                         max_no_correlated = 10)

  delta <- rnorm(J,sd = 2)

  sparsity <- numeric(10)
  maxes <- numeric(10)
  sums <- numeric(10)
  unobs <- numeric(10)
  for(draw in 1:10){
  Y <- generate_sim_data(n_per_group = 25,
                    b0 = c(rep(-8,J/10),seq(-4,-1,length.out = 0.8*J),rep(3,J/10)),
                    b1 = rep(0,J),
                    delta = delta,
                    z_mean = 3,
                    nb_size = .125,
                    cor_structs = cor_structs)
  sparsity[draw] <- mean(Y ==0)
  maxes[draw] <- max(Y)
  sums[draw] <- sum(Y)
  unobs[draw] <- sum(apply(Y, 2, function(x) mean(x ==0))==1)
  }

  expect_equal(mean(sparsity),0.8,tolerance = 0.1)
  expect_true(max(maxes)<1e6 & min(maxes)>.5e4)
  expect_equal(mean(unobs),4.5,tolerance = 1)

  J <- 250
  n_per_group <- 25
  set.seed(0)
  cor_structs <-
    generate_sim_setting(J,
                         n_correlations = 20,
                         max_no_correlated = 10)

  delta <- rnorm(J,sd = 2)
  for(draw in 1:10){
    Y <- generate_sim_data(n_per_group = 25,
                           b0 = c(rep(-8,J/10),seq(-4,-1,length.out = 0.8*J),rep(3,J/10)),
                           b1 = rep(0,J),
                           delta = delta,
                           z_mean = 3,
                           nb_size = .125,
                           cor_structs = cor_structs)
    sparsity[draw] <- mean(Y ==0)
    maxes[draw] <- max(Y)
    sums[draw] <- sum(Y)
    unobs[draw] <- sum(apply(Y, 2, function(x) mean(x ==0))==1)
  }

  expect_equal(mean(sparsity),0.8,tolerance = 0.1)
  expect_true(max(maxes)<5e6 & min(maxes)>1e4)
  expect_equal(mean(unobs)/250,0.07,tolerance = 0.05)


})


test_that("We are actually inducing correlation", {

  J <- 50
  n_per_group <- 50
  set.seed(0)
  cor_structs <-
    generate_sim_setting(J,
                         n_correlations = 10,
                         max_no_correlated = 10)

  delta <- rnorm(J,sd = 2)

  sparsity <- numeric(10)
  maxes <- numeric(10)
  sums <- numeric(10)
  unobs <- numeric(10)
    Y_circ <- generate_sim_data(n_per_group = 500,
                           b0 = c(rep(-9,5),seq(-5,-2,length.out = J - 10),rep(2,5)),
                           b1 = rep(0,J),
                           delta = delta,
                           z_mean = 3,
                           nb_size = .125,
                           cor_structs = cor_structs,
                           return_Y_circ = TRUE)

    expect_equal(min(cor(Y_circ)),-.79,tolerance = 0.05)
    expect_equal(max(cor(Y_circ)[upper.tri(cor(Y_circ))]),0.79,
                 tolerance = 0.05)

})

# test_that("Null simulations for choices of n = 25 and J = 50 work", {
#
#   J <- 60
#   set.seed(0)
#   cor_structs <-
#     generate_sim_setting(J,
#                          n_correlations = 10,
#                          max_no_correlated = 10)
#
#   delta <- rnorm(J,sd = 2)
#
#
#   n25_J60 <-
#     do_one_sim(n_per_group = 25,
#                          J = 60,
#                          delta = delta,
#                          b1_magnitude = 0,
#                          cor_structs = cor_structs,
#                seed = 0,
#                          nboot = 2)
#
#
#   n100_J60 <-
#     do_one_sim(n_per_group = 100,
#                J = 60,
#                delta = delta,
#                b1_magnitude = 0,
#                cor_structs = cor_structs,
#                seed = 0,
#                nboot = 2)
#
#   J <- 240
#   set.seed(0)
#   cor_structs <-
#     generate_sim_setting(J,
#                          n_correlations = 10,
#                          max_no_correlated = 25)
#
#   delta <- rnorm(J,sd = 2)
#
#   n25_J240 <-
#     do_one_sim(n_per_group = 25,
#                J = 240,
#                delta = delta,
#                b1_magnitude = 0,
#                cor_structs = cor_structs,
#                seed = 0,
#                nboot = 2)
#
#   n100_J240 <-
#     do_one_sim(n_per_group = 100,
#                J = 240,
#                delta = delta,
#                b1_magnitude = 0,
#                cor_structs = cor_structs,
#                seed = 0,
#                nboot = 2)
#   expect_true(nrow(n25_J60)== 60*4)
#   expect_true(nrow(n100_J60) == 60*4)
#   expect_true(nrow(n25_J240)== 240*4)
#   expect_true(nrow(n100_J240) == 240*4)
#
# })
