rm(list = ls())

library(mcmc)
library(phylogenetic)
library(microbenchmark)
library(coda)

# Simulating data --------------------------------------------------------------
set.seed(13123)

n <- 5e3
b <- 2
s <- 5

X   <- cbind(rnorm(n))
eps <- rnorm(n, sd = s)
y   <- X*b + eps

# LogLikelihood
loglike <- function(params) {
  sum(log(dnorm((y - X*params[1])/params[2])/params[2] ))
}

# Testing out implementations --------------------------------------------------

# MCMC
ans_phylo <- MCMC(loglike, c(1, 5), nbatch = 100, ub = c(100, 10), lb = c(-100, .5), useCpp = TRUE,
            burnin = 0, scale=.5, thin=1L)

# Metrop
ans_mcmc  <- metrop(loglike, c(1, 5), nbatch = 100, scale=.5)
ans_mcmc  <- mcmc(ans_mcmc$batch)

# Benchmark --------------------------------------------------------------------
microbenchmark(
  # Not using Rcpp
  MCMC_R   = MCMC(
    loglike, c(1, 5), nbatch = 100, ub = c(100, 10), lb = c(-100, .5), useCpp = FALSE,
    burnin = 0, scale=.5, thin=1L),
  # Using Rcpp
  MCMC_cpp = MCMC(
    loglike, c(1, 5), nbatch = 100, ub = c(100, 10), lb = c(-100, .5), useCpp = TRUE,
    burnin = 0, scale=.5, thin=1L),
  # The Metrop algorithm
  metrop   = metrop(loglike, c(1, 5), nbatch = 100, scale=.5),
  unit     = "relative",
  times    = 100L
  )
