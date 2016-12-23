rm(list=ls())

library(phylogenetic)

# Reading data -----------------------------------------------------------------
dag <- read.table("data/pthr11848.dag.txt", sep="\t", header = FALSE,
                  col.names = c("NodeId", "TypeId", "ParentId"))

dat <- read.table("data/pthr11848_sorted.txt", sep = "\t", header = FALSE)
colnames(dat) <- c(sprintf("f%02d",1:(ncol(dat)-1)), "LeafId")

# Input data
Z <- as.matrix(dat[,-4])

# Listing offsprings -----------------------------------------------------------
offsprings <- lapply(dat$LeafId, function(x) {
  x <- dag$NodeId[which(dag$ParentId == x)]
  x <- match(x, dat$LeafId)
  
  x[is.na(x)] <- NULL
  if (length(x)) x
  else integer(0)
})

# Substracting one so we canuse it in C++
offsprings <- lapply(offsprings, function(x) {
  if (length(x)) x-1
  else x
})

# Checking who is parent
noffsprings <- sapply(offsprings, length)

# Parameters -------------------------------------------------------------------

psi     <- c(0.020,0.010)
mu      <- c(0.004,.001)
pi_root <- c(1-0.1,.1)

set.seed(1231)


# Step by step calculation -----------------------------------------------------

S   <- states(ncol(Z))
PSI <- leaf_prob(Z, S, psi, noffsprings)
M   <- gain_loss_prob(mu)
PI  <- root_node_prob(pi_root, S)
Pr <- internal_prob(PSI, M, S, noffsprings, offsprings)

# Doing the same in a single step ---------------------------------------------
ll1   <- LogLike(Z, offsprings, noffsprings, psi, mu, pi_root)

all(ll1$S == S)
all(ll1$PSI == PSI)
all(ll1$M == M)
all(ll1$PI == PI)
ll1$Pr[which(ll1$Pr != Pr, arr.ind = TRUE)[,1],]
Pr[which(ll1$Pr != Pr, arr.ind = TRUE)[,1],]

str(offsprings[which(Pr != ll1$Pr)])