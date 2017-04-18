rm(list = ls())

library(aphylo)

# Simulating data
set.seed(13123)

n <- 2e3
b <- 2
s <- 5

X   <- cbind(rnorm(n))
eps <- rnorm(n, sd = s)
y   <- X*b + eps


loglike <- function(params) {
  sum(log(dnorm((y - X*params[1])/params[2])/params[2] ))
}

ans <- MCMC(loglike, c(1, 5), nbatch = 5e4, ub = c(100, 1e3), lb = c(-100, .5), useCpp = TRUE,
            burnin = 0, scale=.5, thin=100)

summary(ans)
plot(ans)
