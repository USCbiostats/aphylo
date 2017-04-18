# ------------------------------------------------------------------------------
# This script generates data using the DGP of the Annotated Phylogenetic Tree
# model and runs the MLE estimator on each set of data generated, in particular
# N. P sets the number of functions to simulate, n the number of leaf and
# pars the DGP parameters.
# ------------------------------------------------------------------------------

rm(list = ls())

# Simulation parameters
N    <- 1e3
n    <- 1000
P    <- 1
pars <- c(psi0=.05, psi1=.1, mu0=.2, mu1=.06, Pi=.3)

# Preparing the parallel session
library(parallel)
cl <- makeCluster(15)

invisible(clusterEvalQ(cl, {
  library(aphylo)
}))

# Passing parameters to the clusters
clusterExport(cl, c("N", "n", "P", "pars"))

# Doing parLapply
ans <- parSapply(cl, 1:N, function(i) {
  # Simulate a random annotated phylo tree
  tree <- sim_annotated_tree(n, P = P, pars = pars)
  
  # Estimation procedure
  est  <- phylo_mle(rep(.5, 5), tree, useABC = TRUE)

  # Storing the data
  c(est$par, ll=est$ll)
})

stopCluster(cl)

# Checking the results

# > summary(t(ans)[,-6])
#      psi0              psi1              mu0               mu1                Pi        
# Min.   :0.00000   Min.   :0.00000   Min.   :0.02911   Min.   :0.01897   Min.   :0.4917  
# 1st Qu.:0.00000   1st Qu.:0.07068   1st Qu.:0.16726   1st Qu.:0.05346   1st Qu.:0.5000  
# Median :0.05026   Median :0.09726   Median :0.19321   Median :0.06452   Median :0.5000  
# Mean   :0.11410   Mean   :0.14116   Mean   :0.19113   Mean   :0.07435   Mean   :0.5000  
# 3rd Qu.:0.12340   3rd Qu.:0.12166   3rd Qu.:0.22075   3rd Qu.:0.07999   3rd Qu.:0.5000  
# Max.   :1.00000   Max.   :1.00000   Max.   :0.34032   Max.   :0.33260   Max.   :0.5053  
#
# Overall we see that the estimator gets very close to the true parameter: Checking out
# the bias:
#
# > apply(ans, 1, median)[-6] - pars
#         psi0          psi1           mu0           mu1            Pi 
# 0.0002616432 -0.0027352608 -0.0067947610  0.0045169110  0.1999999928 
#
# Around 3 digits precision. All but Pi (which we did suspect returns .5 by construction)
#