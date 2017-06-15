rm(list = ls())

library(aphylo)
library(coda)

# Number of simulations (samples) to draw
nsim <- 13096

# Function to draw parameters:
# This creates 5 + 1 parameters, namely
#  - psi0 Probability of mislabel a 0
#  - psi1 Probability of mislabel a 1
#  - mu0  Probability of gain a 1
#  - mu1  Probability of loss a 1
#  - Pi   Root node probability
#  - % of missings
draw_par <- function() {
  structure(c(rbeta(5, 1, 20), runif(1, .1, .5)), names = c("psi0", "psi1", "mu0", "mu1", "Pi", "missing"))
}

# Function to simulate annotations on a given tree
read_and_sim <- function(fn, p) {
  # Reading
  dat <- ape::read.tree(fn, nmax=1)
  
  # Simulating
  sim_annotated_tree(
    tree = as_po_tree(dat$edge),
    psi = p[1:2],
    mu  = p[3:4],
    Pi  = p[5]
    )
}

# Function to drop annotations
# @param pcent Numeric scalar. Proportion of missings
# @param proportions Numeric vector of length 2. Relative probability that 1s
#  and 0s are kept, i.e., if c(.7, .3) ones will be kept with a 70% chance and
#  0s will be kept with 30% chance.
drop_data <- function(dat, pcent, proportions = c(.2, .8)) {
  
  # Identifying some leaf nodes
  ids <- which(dat$noffspring == 0)
  
  # Which ones we will not drop
  ids <- sample(ids, floor( (1 - pcent)*length(ids)),
                prob = proportions[dat$annotations[ids] + 1]
                )
  
  # Getting the set that will be set to be zero
  ids <- setdiff(1:length(dat$noffspring), ids)
  dat$annotations[ids,] <- 9
  
  dat
}

# Function to estimate model using mle
mle_lite <- function(dat, abc, priors = NULL) {
  # Try to estimate the model
  ans <- if (abc) tryCatch(phylo_mle(dat, method="ABC", priors=priors), error = function(e) e)
    else tryCatch(phylo_mle(dat, priors = priors), error = function(e) e)
  
  
  # If it didn't worked, then return error
  if (inherits(ans, "error")) return(ans)
  
  # ans[c("par", "ll", "counts", "convergence", "message", "method", "varcovar")]
  ans
}

# Function to estimate model using MCMC
mcmc_lite <- function(dat, par, nbatch = 5e5L, nchains = 5L, burnin=1e4, thin=100, priors = NULL) {
  
  # Try to estimate the model
  ans <- tryCatch(phylo_mcmc(
    par, dat,
    control = list(nbatch = nbatch, nchains=nchains, burnin = burnin, thin=thin),
    priors = priors
    ),
    error = function(e) e
    )
  
  # If it didn't worked, then return error
  if (inherits(ans, "error")) return(ans)
  
  # ans[["dat"]] <- NULL
  
  return(ans)
  
  
}

# Function to sample a tree, and simulate annotations on it
tree_files <- list.dirs("../PANTHER11.1/books", full.names = TRUE, recursive = FALSE)
tree_files <- paste0(tree_files, "/tree.tree")

# Checking which ones exists
tree_files_exists <- sapply(tree_files, file.exists)
table(tree_files_exists) # 13096 TRUE

# Sampling files and getting the seed
set.seed(1133)
sample_files <- sample(tree_files, nsim)
parameters   <- lapply(1:nsim, function(x) draw_par())

dat          <- Map(function(fn, p) read_and_sim(fn, p), fn = sample_files, p = parameters)
dat_obs      <- Map(function(d, p) drop_data(d, p[["missing"]]), d = dat, p = parameters)

x <- lapply(dat, "[[", "annotations")
x <- lapply(x, table)
x <- lapply(x, names)
table(unlist(x))

save.image("playground/simulations/data_and_functions.rda")


