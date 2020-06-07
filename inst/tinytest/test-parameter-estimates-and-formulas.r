suppressMessages(library(coda))

# test_that("x ~ mu", {
  
  # Setting the default with no multicore
  pars <- aphylo:::APHYLO_DEFAULT_MCMC_CONTROL
  pars$nsteps  <- 2e3
  pars$burnin  <- 1e3L
  pars$nchains <- 1L
  pars$thin    <- 1L
  pars$conv_checker <- NULL
  pars <- c(pars, list(conv_checker = NULL))
  pars$kernel  <- fmcmc::kernel_am(ub = .9999, lb = 0.001, freq = 1L)

  # Data generating process
  set.seed(7223)
  # x <- rdrop_annotations(raphylo(40), .6)
  x <- raphylo(40)

  mypriors <- function(z) dbeta(z, 2, 10)

  set.seed(1)
  ans0 <- aphylo_mcmc(x ~ mu_d, priors = mypriors, control = pars)
  
  fun <- function(p) {
    ans <- LogLike(
      tree = x,
      psi  = c(0,0),
      mu_d = p[c("mu_d0", "mu_d1")],
      mu_s = p[c("mu_d0", "mu_d1")],
      eta  = -c(1, 1),
      Pi   = p["mu_d0"]/(p["mu_d0"] + p["mu_d1"]),
      verb_ans = FALSE 
      )$ll + sum(log(mypriors(p))) # + + 0.693147180559945 * sum(Nann(x)) * sum(Nannotated(x))

    if (is.infinite(ans))
      ans <- - .Machine$double.xmax * 1e-10
    
    ans
  }

  # Running the raw MCMC
  pars$kernel <- fmcmc::kernel_am(ub = .9999, lb = 0.001, freq = 1L)
  x$reduced_pseq <- x$pseq
  set.seed(1)
  ans1 <- suppressWarnings({
    do.call(
      fmcmc::MCMC, c(
        list(
          fun = fun,
          initial = aphylo:::APHYLO_PARAM_DEFAULT[c("mu_d0", "mu_d1")]
        ),
        pars
      ))
    })
    
  expect_equal(summary(ans1)$statistics[,"Mean"], ans0$par)

# })

# test_that("x ~ mu_d + mu_s + psi + Pi", {
  
  # Data generating process
  set.seed(7223)
  x <- rdrop_annotations(raphylo(40), .6)
  
  mypriors <- function(z) dbeta(z, 2, 10)
  
  pars <- aphylo:::APHYLO_DEFAULT_MCMC_CONTROL
  pars$nsteps  <- 2e3
  pars$burnin  <- 1e3
  pars$thin    <- 10
  pars$nchains <- 1L
  pars$conv_checker <- NULL
  pars <- c(pars, list(conv_checker = NULL))
  pars$kernel  <- fmcmc::kernel_am(ub = .9999, lb = 0.001, freq = 1L)
  
  set.seed(1)
  dflts <- aphylo:::APHYLO_PARAM_DEFAULT[-c(7:8)]
  ans0 <- suppressWarnings(
    aphylo_mcmc(x ~ mu_d + mu_s + psi + Pi, params = dflts, priors = mypriors, control = pars)
  )
  
  mcmc_0 <- as.list(fmcmc::LAST_MCMC)
  
  fun <- function(p) {
    ans <- aphylo::LogLike(
      tree = x,
      psi  = p[c("psi0", "psi1")],
      mu_d = p[c("mu_d0", "mu_d1")],
      mu_s = p[c("mu_s0", "mu_s1")],
      eta  = -c(1, 1),
      Pi   = p["Pi"],
      verb_ans = FALSE 
    )$ll + sum(log(mypriors(p))) # log(2^prod(dim(x$tip.annotation)))
    
    if (!is.finite(ans))
      ans <- -.Machine$double.xmax * 1e-10
    
    ans
  }
  
  # Running the raw MCMC
  pars <- aphylo:::APHYLO_DEFAULT_MCMC_CONTROL
  pars$nsteps  <- 2e3
  pars$burnin  <- 1e3
  pars$nchains <- 1L
  pars$thin    <- 10
  pars$conv_checker <- NULL
  pars <- c(pars, list(conv_checker = NULL))
  set.seed(1)
  ans1 <- suppressWarnings({
    
    do.call(
      fmcmc::MCMC, c(
        list(
          fun = fun,
          initial = aphylo:::APHYLO_PARAM_DEFAULT[-c(7:8)],
          kernel  = fmcmc::kernel_am(ub = .9999, lb = 0.001, freq = 1L)
          ),
        pars
      ))
    
    })
  
  mcmc_1 <- as.list(fmcmc::LAST_MCMC)
  
  expect_equal(summary(ans1)$statistics[,"Mean"], ans0$par)
  # ans0$fun(colMeans(ans1), dat = ans0$dat, priors = ans0$priors)
  # fun(colMeans(ans1))
# })

