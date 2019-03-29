context("The formulas are writing the right model")

suppressMessages(library(coda))

test_that("x ~ mu", {
  
  # Setting the default with no multicore
  pars <- aphylo:::APHYLO_DEFAULT_MCMC_CONTROL
  pars$multicore <- FALSE

  # Data generating process
  set.seed(7223)
  x <- rdrop_annotations(raphylo(40), .6)

  mypriors <- function(z) dbeta(z, 2, 10)

  set.seed(1)
  ans0 <- suppressWarnings(aphylo_mcmc(x ~ mu, priors = mypriors, control = pars))
  

  fun <- function(p) {
    ans <- LogLike(
      tree = x,
      psi  = c(0,0),
      mu   = p[c("mu0", "mu1")],
      eta  = c(.5, .5),
      Pi   = p["mu0"]/(p["mu0"] + p["mu1"]),
      verb_ans = FALSE 
      )$ll + sum(log(mypriors(p))) + log(2^prod(dim(x$tip.annotation)))

    if (is.infinite(ans))
      ans <- .Machine$double.xmax * sign(ans) * 1e-10
    
    ans
  }

  # Running the raw MCMC
  set.seed(1)
  
  ans1 <- suppressWarnings({
    do.call(
      amcmc::MCMC, c(
        list(
          fun = fun,
          initial = aphylo:::APHYLO_PARAM_DEFAULT[c("mu0", "mu1")]
        ),
        pars
      ))
    })
    
  expect_equal(summary(ans1)$statistics[,"Mean"], ans0$par)

})

test_that("x ~ mu + psi + Pi", {
  
  # Data generating process
  set.seed(7223)
  x <- rdrop_annotations(raphylo(40), .6)
  
  mypriors <- function(z) dbeta(z, 2, 10)
  
  set.seed(1)
  ans0 <- suppressWarnings(aphylo_mcmc(x ~ mu + psi + Pi, priors = mypriors,
                                       control = list(multicore=FALSE)))
  
  
  fun <- function(p) {
    ans <- aphylo::LogLike(
      tree = x,
      psi  = p[c("psi0", "psi1")],
      mu   = p[c("mu0", "mu1")],
      eta  = c(.5, .5),
      Pi   = p["Pi"],
      verb_ans = FALSE 
    )$ll + sum(log(mypriors(p))) + log(2^prod(dim(x$tip.annotation)))
    
    if (is.infinite(ans))
      ans <- .Machine$double.xmax * sign(ans) * 1e-10
    
    ans
  }
  
  # Running the raw MCMC
  pars <- aphylo:::APHYLO_DEFAULT_MCMC_CONTROL
  pars$multicore <- FALSE
  # pars$nchains   <- 1L
  set.seed(1)
  ans1 <- suppressWarnings({
    
    do.call(
      amcmc::MCMC, c(
        list(
          fun = fun,
          initial = aphylo:::APHYLO_PARAM_DEFAULT[-c(5:6)]
          ),
        pars
      ))
    
    })
  
  expect_equal(summary(ans1)$statistics[,"Mean"], ans0$par)
  
})