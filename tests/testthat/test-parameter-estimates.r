context("Parameter estimation")

set.seed(154)
n <- 100
P <- 1
psi <- c(.01, .01)
mu  <- c(.02, .05)
Pi  <- .2

dat <- raphylo(n, P=P, psi = psi, mu = mu, Pi = Pi)

# Estimation via L-BFGS-B

ans0 <- suppressWarnings({
  aphylo_mle(
    dat ~ mu + psi + eta + Pi,
    params = c(.05, .05, .05, .05, .9, .9, .5), lower = 0, upper=1)
})

# Methods ----------------------------------------------------------------------
test_that("Methods", {
  
  # Printing
  expect_output(suppressWarnings(print(ans0)), "ESTIMATION OF")
  
  # Plotting
  expect_silent(plot(ans0))
  expect_silent(plot_logLik(ans0))
  
  # Extracting coef and others
  expect_equal(coef(ans0), ans0$par)
  expect_equal(vcov(ans0), ans0$varcovar)
  
  
})

# ------------------------------------------------------------------------------
test_that("MCMC", {
  ans1 <- suppressWarnings(
    aphylo_mcmc(dat ~ mu + psi + eta + Pi, params = ans0$par,
                priors = function(p) 1,
                control = list(nsteps = 1e4, burnin=5e3, thin=20)
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
})

test_that("MCMC: in a degenerate case all parameters goes to the prior", {
  
  
  set.seed(1)
  dat <- suppressWarnings(raphylo(50, Pi=0, mu=c(0, 0), psi=c(0,0)))
  dat$tip.annotation[] <- 9L
  
  ans1 <- suppressWarnings(
    aphylo_mcmc(dat ~ mu + psi + eta(0,1) + Pi,
                params = c(rep(2/12, 4), .5, .5,2/12),
                priors = function(x) dbeta(x, 2, 10),
                control = list(
                  nsteps = 4e4, burnin=1e4, nchains=2,
                  kernel = fmcmc::kernel_reflective(lb = 0, ub = 1, scale = .05)
                  ),
                check_informative = FALSE
                )
    )
  
  ans2 <- suppressWarnings(
    aphylo_mcmc(
      dat ~ mu + psi + eta(0,1) + Pi,
      params = c(rep(2/22, 4), .5,.5,2/22), 
      priors = function(x) dbeta(x, 2, 20),
      control = list(
        nsteps = 4e4, burnin=1e4, nchains=2,
        kernel = fmcmc::kernel_reflective(lb = 0, ub = 1, scale = .05)
        ),
      check_informative = FALSE
      )
  )
  
  
  # Should converge to the prior
  expect_equal(unname(coef(ans1))[-c(5:6)], rep(2/12, 5), tol=.025)
  expect_equal(unname(coef(ans2))[-c(5:6)], rep(2/22, 5), tol=.025)
  
})
