rm(list = ls())

library(mcmc)
library(MCMCpack)
library(phylogenetic)

data(experiment)
data(tree)

# Preprocessing the data
dat <- get_offspring(experiment, "LeafId",
                     tree, "NodeId", "ParentId")

priors <- function(params) {
  dbeta(params[1:2], 1, 9)
}

expit <- function(x)
  exp(x) / (1 + exp(x))

fun <- function(params) {
  if (any(params > 1 | params < 0))
    return(-Inf)
  
  psi <- params[1:2]
  mu  <- params[3:4]
  Pi  <- params[5]
  Pi  <- c(1 - Pi, Pi)
  
  LogLike(dat$experiment,
          dat$offspring,
          dat$noffspring,
          psi,
          mu,
          Pi,
          FALSE)$ll +
    sum(log(priors(params)))
}

# Metropolis Random Walk -------------------------------------------------------

# Running the algorithm
set.seed(1231)
ans <- metrop(fun, rep(.5, 5), nbatch = 2e3, scale = 0.25)

# Checking answer
sprintf("%.5f" , ans$initial)
sprintf("%.5f" , ans$final)
fun(ans$final)
ans$accept

boxplot(ans$batch)

# Diagnostics
tail(ans$batch)


plot(ts(ans$batch))
acf(ans$batch)

# Tempering algorithm ----------------------------------------------------------

fun <- function(params) {
  # Getting the i
  i      <- params[1]
  params <- params[-1]
  
  if (any(params > 1 | params < 0))
    return(-Inf)
  
  psi <- params[1:2]
  mu  <- params[3:4]
  Pi  <- params[5]
  Pi  <- c(1 - Pi, Pi)
  
  - (
    -LogLike(
      dat$experiment,
      dat$offspring,
      dat$noffspring,
      psi,
      mu,
      Pi,
      FALSE
    )$ll -
      sum(log(priors(params)))
  ) ^ (1 / i)
}

ncomp <- 5
neighbors <- matrix(FALSE, ncomp, ncomp)
neighbors[row(neighbors) == col(neighbors) + 1] <- TRUE
neighbors[row(neighbors) == col(neighbors) - 1] <- TRUE

thetas <- c(1L, rep(0.5, 5))

set.seed(123)
ans <- temper(fun,
              thetas,
              neighbors = neighbors,
              nbatch = 1e4,
              scale = .25,
              blen = 10)

boxplot(ans$batch)
plot(ts(ans$batch))
summary(ans$batch)
ans$accepti
ans$acceptx

fun(c(1L,ans$final))
# Metropolist ------------------------------------------------------------------

fun <- function(params) {
  params <- expit(params)
  
  # if (any(params > 1 | params < 0)) return(-1e5)
  
  psi <- params[1:2]
  mu  <- params[3:4]
  Pi  <- params[5]
  Pi  <- c(1 - Pi, Pi)
  
  LogLike(dat$experiment,
          dat$offspring,
          dat$noffspring,
          psi,
          mu,
          Pi,
          FALSE)$ll +
    sum(log(priors(params)))
}


# This method won't run because non-negative definite matrix
set.seed(1231)
ans <- MCMCmetrop1R(
  fun,
  rep(.5, 5),
  optim.lower = 1e-5,
  optim.upper = 1 - 1e-5,
  optim.method = "L-BFGS-B"
)
