# context("Parameter estimation")

set.seed(15)
n <- 200
P <- 1
psi <- c(.01, .01)
mu_d  <- c(.05, .05)
Pi  <- .2

dat <- raphylo(n, P=P, psi = psi, mu_d = mu_d, Pi = Pi)

# Estimation via L-BFGS-B

ans0 <- suppressWarnings({
  aphylo_mle(
    dat ~ mu_d + psi + eta + Pi,
    params = c(.05, .05, .05, .05, .9, .9, .5), lower = 0, upper=1)
})

# Methods ----------------------------------------------------------------------
# test_that("Methods", {
  
  # Printing
  expect_silent(print(ans0))
  
  # Plotting
  expect_silent(plot(ans0))
  expect_silent(plot_logLik(ans0))
  
  # Extracting coef and others
  expect_equal(coef(ans0), ans0$par)
  expect_equal(vcov(ans0), ans0$varcovar)
  
  
# })

# ------------------------------------------------------------------------------
# test_that("MCMC", {
  ans1 <- suppressWarnings(
    aphylo_mcmc(dat ~ mu_d + psi + eta + Pi, params = ans0$par,
                priors = function(p) 1,
                control = list(
                  nsteps = 1e4, burnin=5e3, thin=20,
                  conv_checker = NULL, nchains=1
                  )
                )
  )
  
  # Checking expectations. It is difficult to do so with PI, so we exclude it.
  expect_equal(ans0$par[-7], ans1$par[-7], tolerance = 0.1)
  
  # Plotting
  expect_silent(plot(ans1))
  expect_silent(plot_logLik(ans1))
  
  # Extracting coef and others
  expect_equal(coef(ans1), ans1$par)
  expect_equal(vcov(ans1), ans1$varcovar)
# })

# test_that("MCMC: in a degenerate case all parameters goes to the prior", {
  
  
set.seed(1)
dat <- suppressWarnings(raphylo(50, Pi=0, mu_d=c(0, 0), psi=c(0,0)))
dat$tip.annotation[] <- 9L

ans1 <- suppressWarnings(
  aphylo_mcmc(dat ~ mu_d + psi + eta(0,1) + Pi,
              params = c(rep(2/12, 4), .5, .5,2/12),
              priors = function(x) dbeta(x, 2, 10),
              control = list(
                nsteps = 4e4, burnin=1e4, nchains=2,
                kernel = fmcmc::kernel_normal_reflective(
                  lb = 0, ub = 1, scale = .05,
                  order = "random"
                  ),
                conv_checker = NULL
                ),
              check_informative = FALSE
              )
  )

ans2 <- suppressWarnings(
  aphylo_mcmc(
    dat ~ mu_d + psi + eta(0,1) + Pi,
    params = c(rep(2/22, 4), .5,.5,2/22), 
    priors = function(x) dbeta(x, 2, 20),
    control = list(
      nsteps = 4e4, burnin=1e4, nchains=2,
      kernel = fmcmc::kernel_normal_reflective(
        lb = 0, ub = 1, scale = .05,
        order = "random"),
      conv_checker = NULL
      ),
    check_informative = FALSE
    )
)
  
  
  # Should converge to the prior
expect_true(all(abs(unname(coef(ans1))[-c(5:6)] - rep(2/12, 5)) < .025))
expect_true(all(abs(unname(coef(ans2))[-c(5:6)] - rep(2/22, 5)) < .025))
  
# })

# test_that("mu_d and mu_s work", {
  
  dat   <- c(1, 2, 1, 3, 2, 4, 2, 5, 3, 6, 6, 7, 6, 8, 4, 9, 4, 10)
  dat   <- as.phylo(matrix(dat, ncol=2, byrow = TRUE))
  types <- c(1, 1, 1, 1, 1, 0, 1, 1, 1, 1)

  set.seed(1)
  
  tree <- raphylo(
    tree      = dat,
    tip.type  = types[dat$tip.label],
    node.type = types[dat$node.label],
    psi       = c(0, 0),
    eta       = c(1,1),
    mu_d      = c(.9,0),
    mu_s      = c(0,0),
    Pi        = 0
  )
  
  plot(tree)
  
  ans_mcmc <- suppressWarnings(aphylo_mcmc(tree ~ mu_d + mu_s, priors = function(p) {
    c(dbeta(p[c("mu_d0", "mu_d1")], 10, 2), 
      dbeta(p[c("mu_s0", "mu_s1")], 2, 10))
  }))
  
  ans_mle <- aphylo_mle(tree ~ mu_d + mu_s, priors = function(p) {
    c(dbeta(p[c("mu_d0", "mu_d1")], 10, 2), 
      dbeta(p[c("mu_s0", "mu_s1")], 2, 10))
  })
  
  expect_equivalent(predict(ans_mle), predict(ans_mcmc), tol=.15)
  
  
# })
