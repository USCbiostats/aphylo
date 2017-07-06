set.seed(18888)

load("playground/simulations/mcmc_right_prior_estimates.rda")

N   <- 1e3
mcmc_right_prior_estimates_sample <- ans_MCMC_right_prior[1:N]

save(mcmc_right_prior_estimates_sample,
    file = "playground/simulations/mcmc_right_prior_estimates_sample.rda")
