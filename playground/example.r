rm(list=ls())

library(phylogenetic)

# Reading data -----------------------------------------------------------------
data("experiment")
data("tree")

# Listing offspring -----------------------------------------------------------
O <- get_offspring(
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
PSI <- leaf_prob(Z, S, psi, O$noffspring)
M   <- gain_loss_prob(mu)
PI  <- root_node_prob(pi_root, S)
Pr  <- internal_prob(PSI, M, S, O$noffspring, O$offspring)

# Doing the same in a single step ---------------------------------------------
ll1   <- LogLike(Z, O$offspring, O$noffspring, psi, mu, pi_root)

all(ll1$S == S)
all(ll1$PSI == PSI)
all(ll1$M == M)
all(ll1$PI == PI)
ll1$Pr[which(ll1$Pr != Pr, arr.ind = TRUE)[,1],]
Pr[which(ll1$Pr != Pr, arr.ind = TRUE)[,1],]

str(O$offspring[which(Pr != ll1$Pr)])

# Optimization -----------------------------------------------------------------

fun <- function(params) {
  psi     <- params[1:2] # c(0.020,0.010)
  mu      <- params[3:4] # c(0.004,.001)
  pi_root <- c(params[5], 1-params[5]) # c(1-0.1,.1) 
  
  -LogLike(Z, O$offspring, O$noffspring, psi, mu, pi_root)$ll
}

# A (not that) bad example
set.seed(123)
library(ABCoptim)
seedpar <- runif(5)
ans_abcoptim <- abc_cpp(seedpar, fun, 1, 0, maxCycle = 500, criter = 50)
message(
  "ABCoptim results",
  "\n - Params : ",
  paste0(sprintf("%06.4f", ans_abcoptim$par), collapse=", "),
  "\n - ll     : ", sprintf("%f", -ans_abcoptim$value)
)


# Using the -numDerive- package
library(numDeriv)

expit <- function(x) exp(x)/(1 + exp(x))
logit <- function(x) log(x/(1-x))

fun <- function(params) {
  params  <- expit(params)
  psi     <- params[1:2] # c(0.020,0.010)
  mu      <- params[3:4] # c(0.004,.001)
  pi_root <- c(params[5], 1-params[5]) # c(1-0.1,.1) 

  -LogLike(Z, O$offspring, O$noffspring, psi, mu, pi_root)$ll
}

niter <- 20
PARAMS <- matrix(ncol=5, nrow=niter)
params <- logit(ans_abcoptim$par)
criter <- 1e-5
for (i in 1:niter) {
  
  PARAMS[i,] <- params
  params0    <- params
  
  # Computing jacobian and hessian
  fun_jacb   <- jacobian(fun, params, method.args=list(d=.05))
  fun_hess   <- hessian(fun, params, method.args=list(d=.05))
  
  # Updating step
  params <- params - fun_jacb %*% solve(fun_hess, tol = 1e-20)
  
  # Error
  if (is.na(fun(params))) 
    stop("Undefined value of fun(params).")
  
  # Stopping criteria
  val <- abs(fun(params) - fun(PARAMS[i,]))
  if (val < criter) {
    message(
      "NR Results",
      "\n - Params : ",
      paste0(sprintf("%06.4f",expit(params)), collapse=", "),
      "\n - ll     : ", sprintf("%f", -fun(params))
      )
    break
  }

  print(fun(params))
  
}

plot.new()
plot.window(ylim=c(0,1),xlim=c(1,i))
for (j in 1:i)
  lines(data.frame(1:i,PARAMS[1:i,j]), col=j)
axis(2);axis(1)

