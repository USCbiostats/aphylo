rm(list=ls())

library(phylogenetic)

# Reading data -----------------------------------------------------------------
data("experiment")
data("tree")

# Listing offsprings -----------------------------------------------------------
O <- get_offsprings(
  experiment, "LeafId", 
  tree, "NodeId", "ParentId"
)

# Input data
Z <- as.matrix(O$experiment[,-4])

# Parameters -------------------------------------------------------------------
psi     <- c(0.020,0.010)
mu      <- c(0.004,.001)
pi_root <- c(1-0.1,.1)

set.seed(1231)


# Step by step calculation -----------------------------------------------------

S   <- states(ncol(Z))
PSI <- leaf_prob(Z, S, psi, O$noffsprings)
M   <- gain_loss_prob(mu)
PI  <- root_node_prob(pi_root, S)
Pr  <- internal_prob(PSI, M, S, O$noffsprings, O$offsprings)

# Doing the same in a single step ---------------------------------------------
ll1   <- LogLike(Z, O$offsprings, O$noffsprings, psi, mu, pi_root)

all(ll1$S == S)
all(ll1$PSI == PSI)
all(ll1$M == M)
all(ll1$PI == PI)
ll1$Pr[which(ll1$Pr != Pr, arr.ind = TRUE)[,1],]
Pr[which(ll1$Pr != Pr, arr.ind = TRUE)[,1],]

str(O$offsprings[which(Pr != ll1$Pr)])

# Optimization -----------------------------------------------------------------

fun <- function(params) {
  psi     <- params[1:2] # c(0.020,0.010)
  mu      <- params[3:4] # c(0.004,.001)
  pi_root <- c(params[5], 1-params[5]) # c(1-0.1,.1) 
  
  -LogLike(Z, O$offsprings, O$noffsprings, psi, mu, pi_root)$ll
}

# A bad example
set.seed(123)
library(ABCoptim)
ans <- abc_cpp(runif(5), fun, 1, 0, maxCycle = 500, criter = 100);tail(ans$hist)

# Using the -numDerive- package
library(numDeriv)

trunfun <- function(x) {
  ifelse(x>1, 1, ifelse(x<0, 0, x))
}

fun <- function(params) {
  psi     <- params[1:2] # c(0.020,0.010)
  mu      <- params[3:4] # c(0.004,.001)
  pi_root <- c(params[5], 1-params[5]) # c(1-0.1,.1) 
  lambdas <- params[6:10]

  -LogLike(Z, O$offsprings, O$noffsprings, psi, mu, pi_root)$ll - sum(lambdas*params[-c(6:10)])
}

set.seed(123)
niter <- 20
PARAMS <- matrix(ncol=10, nrow=niter)
params <- runif(10)
for (i in 1:niter) {
  
  PARAMS[i,] <- params
  params0    <- params
  fun_jacb   <- jacobian(fun, params, method.args=list(d=.005))
  fun_hess   <- hessian(fun, params, method.args=list(d=.005))
  if (is.na(fun(params))) break
  params <- params - fun_jacb %*% solve(fun_hess, tol = 1e-50)

  print(fun(params))
  
}

plot.new()
plot.window(ylim=c(0,1),xlim=c(1,i))
for (j in 1:i)
  lines(data.frame(1:i,PARAMS[1:i,j]), col=j)
axis(2);axis(1)

