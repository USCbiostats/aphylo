rm(list =ls())


# In this example we estimate the parameter for a dataset with ----------------
# With 5,000 draws from a MVN() with parameters M and S.

# Loading the required packages
library(mvtnorm)
library(aphylo)
library(coda)

# Parameters and data simulation
S <- cbind(c(.8, .2), c(.2, 1))
M <- c(0, 1)

set.seed(123)
D <- rmvnorm(5e3, mean = M, sigma = S)

# Function to pass to MCMC
fun <- function(pars) {
  # Putting the parameters in a sensible way
  m <- pars[1:2]
  s <- cbind( c(pars[3], pars[4]), c(pars[4], pars[5]) )
  
  # Computing the unnormalized log likelihood
  sum(log(dmvnorm(D, m, s)))
}

# Calling MCMC
ans <- MCMC(
  fun,
  initial = c(mu0=5, mu1=5, s0=5, s01=0, s2=5), 
  lb      = c(-10, -10, .01, -5, .01),
  ub      = 5,
  nbatch  = 1e5,
  thin    = 20,
  scale   = .01,
  burnin  = 5e3,
  useCpp  = TRUE
)

# Checking out the outcomes
plot(ans)
summary(ans)


