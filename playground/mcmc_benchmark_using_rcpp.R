rm(list =ls())

library(phylogenetic)
library(microbenchmark)

set.seed(1231)

n <- 1e3
pars <- c(mean = 2.6, sd = 3)

# Generating data and writing the log likelihood function
D <- rnorm(n, pars[1], pars[2])
fun <- function(x) {
  x <- log(dnorm(D, x[1], x[2]))
  if (any(is.infinite(x)))
    return(-Inf)
  sum(x)
}

set.seed(1231)
ans_R <- MCMC(fun, c(mu=1, sigma=1), nbatch = 2e3, scale = .1, ub = 10, lb = 0)
set.seed(1231)
ans_cpp <- MCMC(fun, c(mu=1, sigma=1), nbatch = 2e3, scale = .1, ub = 10, lb = 0, useCpp = TRUE)

summary(ans_R)
summary(ans_cpp)

microbenchmark(
  ans_R   = MCMC(fun, c(mu=1, sigma=1), nbatch = 5e2, burnin = 0, scale = .1, ub = 10, lb = 0),
  ans_cpp = MCMC(fun, c(mu=1, sigma=1), nbatch = 5e2, burnin = 0, scale = .1, ub = 10, lb = 0, useCpp = TRUE),
  unit = "relative", times = 1000
)

# Unit: relative
#    expr      min       lq     mean   median       uq      max neval
#   ans_R 1.461158 1.461158 1.461158 1.461158 1.461158 1.461158     1
# ans_cpp 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000     1
