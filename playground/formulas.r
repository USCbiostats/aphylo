library(aphylo)

eta <- function(eta0 = NULL, eta1 = NULL, equal=FALSE) {
  
  
  
}

# Include all parameters
~ psi + mu + eta + Pi

# Pi is a function of mu
~ psi + mu + eta

# There's no publication bias
~ psi + mu + Pi

# No publication bias, and Pi is fixed to .5
~ psi + mu + Pi(TRUE)

# No publication bias, and Pi is fixed to .5
# and psi0 is fixed to .1
~ psi(fixed = c(0, 1)) + mu + Pi(fixed = TRUE)

list(
  
  fixed = structure(rep(FALSE, 7), names = c("psi0", "psi1", "mu0", "mu1", "eta0", "eta1", "Pi"))
  
)