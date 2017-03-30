context("Parameter estimation")

set.seed(121)
n <- 50
P <- 2
psi <- c(.01, .01)
mu  <- c(.04, .02)
Pi  <- .5

dat <- sim_annotated_tree(n, P=P, psi = psi, mu = mu, Pi = Pi)

# Estimation via L-BFGS-B
ans0 <- phylo_mle(dat)

# ------------------------------------------------------------------------------
test_that("MLE using ABC or L-BFGS-B are close", {
  
  # Estimating models
  ans1 <- phylo_mle(dat, method="ABC")
  
  # Checking expectations
  expect_equal(ans0$par, ans1$par, tolerance = 0.025)
})

# ------------------------------------------------------------------------------
test_that("MCMC", {
  ans2 <- suppressWarnings(
    phylo_mcmc(params = ans0$par, dat, control = list(nbatch = 1e4, burnin=500, thin=20))
  )
  
  # Checking expectations
  expect_equal(ans0$par, ans2$par, tolerance = 0.1)
})

