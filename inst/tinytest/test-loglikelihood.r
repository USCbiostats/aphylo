# context("LogLikelihood function")

data("faketree")
data("fakeexperiment")

# Adding an NA
# fakeexperiment[1,2] <- NA

# Parameters
psi    <- c(.01, .02)
mu     <- c(.05, .1)
eta    <- c(.9, .7)
Pi     <- 1 - .3
errtol <- 1e-15

O  <- new_aphylo(tip.annotation = fakeexperiment[,-1L], tree = as.phylo(faketree))

S  <- states(2)
Pr <- LogLike(
  new_aphylo_pruner(O), psi = psi, mu_d = mu, mu_s = mu, eta=eta,
  Pi = Pi)$Pr[[1]]
  
  # aphylo:::probabilities(
  # with(O, rbind(tip.annotation, node.annotation)),
  # O$pseq, psi, mu, eta, S, O$offspring)

# Pr1 <- LogLike(O, psi = psi, mu =mu,eta=eta, Pi = Pi)
# Pr2 <- LogLike(new_aphylo_pruner(O), psi = psi, mu =mu,eta=eta, Pi = Pi)
# 
# Pr2$Pr - Pr1$Pr[[1]]


# Checking Leaf Probabilities --------------------------------------------------
# test_that("Leaf Probabilities", {
  # Checking cases compared to the states (0,0) (1,0) (0,1) (1,1)
  PrRaw   <- Pr
  PrRaw[] <- 1
  
  # Node 3 (0,0)
  PrRaw[4, 1] <- (1 - psi[1])^2 * eta[1]^2
  PrRaw[4, 2] <- psi[2] * (1 - psi[1]) * eta[1]^2
  PrRaw[4, 3] <- (1 - psi[1]) * psi[2] * eta[1]^2
  PrRaw[4, 4] <- psi[2]^2 * eta[1]^2
  
  # Node 4 (0,1)
  PrRaw[5, 1] <- (1 - psi[1])*psi[1] * eta[1] * eta[2]
  PrRaw[5, 2] <- psi[2] * psi[1] * eta[1] * eta[2]
  PrRaw[5, 3] <- (1 - psi[1]) * (1 - psi[2]) * eta[1] * eta[2]
  PrRaw[5, 4] <- psi[2] * (1 - psi[2]) * eta[1] * eta[2]
  
  # Node 5 (1,0)
  PrRaw[6, 1] <- psi[1]*(1 - psi[1]) * eta[2] * eta[1] 
  PrRaw[6, 2] <- (1 - psi[2]) * (1 - psi[1]) * eta[2] * eta[1] 
  PrRaw[6, 3] <- psi[1] * psi[2] * eta[2] * eta[1] 
  PrRaw[6, 4] <- (1 - psi[2]) * psi[2] * eta[2] * eta[1] 
  
  # Node 6 (1,1)
  PrRaw[7, 1] <- psi[1] ^ 2 * eta[2]^2
  PrRaw[7, 2] <- (1 - psi[2]) * psi[1] * eta[2]^2
  PrRaw[7, 3] <- psi[1] * (1 - psi[2]) * eta[2]^2
  PrRaw[7, 4] <- (1 - psi[2]) ^ 2 * eta[2]^2
  
  # These should be identical
  expect_equivalent(PrRaw[4:7,], Pr[c(5:7, 1:4),][4:7,])
# })


# Checking Internal Probabilities ----------------------------------------------
# test_that("Internal Probabilities", {
  M  <- aphylo:::prob_mat(mu)
  # Checking cases compared to the states (0,0) (1,0) (0,1) (1,1)
  
  PrRaw <- Pr[c(5:7, 1:4),]
  PrRaw[1:3,] <- 1
  
  # Node 1 (0,0)
  of1 <- 
    PrRaw[4:5,1]*M[1,1]*M[1,1] + # (0,0), (0,0)
    PrRaw[4:5,2]*M[1,2]*M[1,1] + # (0,0), (1,0)
    PrRaw[4:5,3]*M[1,1]*M[1,2] + # (0,0), (0,1)
    PrRaw[4:5,4]*M[1,2]*M[1,2]   # (0,0), (1,1)
  
  PrRaw[2,1] <- prod(of1)
  
  # Node 2 (0,0)
  of1 <- 
    PrRaw[6:7,1]*M[1,1]*M[1,1] + # (0,0), (0,0)
    PrRaw[6:7,2]*M[1,2]*M[1,1] + # (0,0), (1,0)
    PrRaw[6:7,3]*M[1,1]*M[1,2] + # (0,0), (0,1)
    PrRaw[6:7,4]*M[1,2]*M[1,2]   # (0,0), (1,1)
  
  PrRaw[3,1] <- prod(of1)
  
  # Node 1 (1,0)
  of1 <- 
    PrRaw[4:5,1]*M[2,1]*M[1,1] + # (1,0), (0,0)
    PrRaw[4:5,2]*M[2,2]*M[1,1] + # (1,0), (1,0)
    PrRaw[4:5,3]*M[2,1]*M[1,2] + # (1,0), (0,1)
    PrRaw[4:5,4]*M[2,2]*M[1,2]   # (1,0), (1,1)
  
  PrRaw[2,2] <- prod(of1)
  
  # Node 2 (1,0)
  of1 <- 
    PrRaw[6:7,1]*M[2,1]*M[1,1] + # (1,0), (0,0)
    PrRaw[6:7,2]*M[2,2]*M[1,1] + # (1,0), (1,0)
    PrRaw[6:7,3]*M[2,1]*M[1,2] + # (1,0), (0,1)
    PrRaw[6:7,4]*M[2,2]*M[1,2]   # (1,0), (1,1)
  
  PrRaw[3,2] <- prod(of1)
  
  # Node 1 (0,1)
  of1 <- 
    PrRaw[4:5,1]*M[1,1]*M[2,1] + # (0,1), (0,0)
    PrRaw[4:5,2]*M[1,2]*M[2,1] + # (0,1), (1,0)
    PrRaw[4:5,3]*M[1,1]*M[2,2] + # (0,1), (0,1)
    PrRaw[4:5,4]*M[1,2]*M[2,2]   # (0,1), (1,1)
  
  PrRaw[2,3] <- prod(of1)
  
  # Node 2 (0,1)
  of1 <- 
    PrRaw[6:7,1]*M[1,1]*M[2,1] + # (0,1), (0,0)
    PrRaw[6:7,2]*M[1,2]*M[2,1] + # (0,1), (1,0)
    PrRaw[6:7,3]*M[1,1]*M[2,2] + # (0,1), (0,1)
    PrRaw[6:7,4]*M[1,2]*M[2,2]   # (0,1), (1,1)
  
  PrRaw[3,3] <- prod(of1)
  
  # Node 1 (1,1)
  of1 <- 
    PrRaw[4:5,1]*M[2,1]*M[2,1] + # (1,1), (0,0)
    PrRaw[4:5,2]*M[2,2]*M[2,1] + # (1,1), (1,0)
    PrRaw[4:5,3]*M[2,1]*M[2,2] + # (1,1), (0,1)
    PrRaw[4:5,4]*M[2,2]*M[2,2]   # (1,1), (1,1)
  
  PrRaw[2,4] <- prod(of1)
  
  # Node 2 (1,1)
  of1 <- 
    PrRaw[6:7,1]*M[2,1]*M[2,1] + # (1,1), (0,0)
    PrRaw[6:7,2]*M[2,2]*M[2,1] + # (1,1), (1,0)
    PrRaw[6:7,3]*M[2,1]*M[2,2] + # (1,1), (0,1)
    PrRaw[6:7,4]*M[2,2]*M[2,2]   # (1,1), (1,1)
  
  PrRaw[3,4] <- prod(of1)
  
  # Node 0 (0,0)
  of1 <- 
    PrRaw[2:3,1]*M[1,1]*M[1,1] + # (0,0), (0,0)
    PrRaw[2:3,2]*M[1,2]*M[1,1] + # (0,0), (1,0)
    PrRaw[2:3,3]*M[1,1]*M[1,2] + # (0,0), (0,1)
    PrRaw[2:3,4]*M[1,2]*M[1,2]   # (0,0), (1,1)
  
  PrRaw[1,1] <- prod(of1)
  
  # Node 0 (1,0)
  of1 <- 
    PrRaw[2:3,1]*M[2,1]*M[1,1] + # (1,0), (0,0)
    PrRaw[2:3,2]*M[2,2]*M[1,1] + # (1,0), (1,0)
    PrRaw[2:3,3]*M[2,1]*M[1,2] + # (1,0), (0,1)
    PrRaw[2:3,4]*M[2,2]*M[1,2]   # (1,0), (1,1)
  
  PrRaw[1,2] <- prod(of1)
  
  # Node 0 (0,1)
  of1 <- 
    PrRaw[2:3,1]*M[1,1]*M[2,1] + # (0,1), (0,0)
    PrRaw[2:3,2]*M[1,2]*M[2,1] + # (0,1), (1,0)
    PrRaw[2:3,3]*M[1,1]*M[2,2] + # (0,1), (0,1)
    PrRaw[2:3,4]*M[1,2]*M[2,2]   # (0,1), (1,1)
  
  PrRaw[1,3] <- prod(of1)
  
  # Node 0 (1,1)
  of1 <- 
    PrRaw[2:3,1]*M[2,1]*M[2,1] + # (1,1), (0,0)
    PrRaw[2:3,2]*M[2,2]*M[2,1] + # (1,1), (1,0)
    PrRaw[2:3,3]*M[2,1]*M[2,2] + # (1,1), (0,1)
    PrRaw[2:3,4]*M[2,2]*M[2,2]   # (1,1), (1,1)
  
  PrRaw[1,4] <- prod(of1)
  
  expect_equal(Pr[c(5:7, 1:4),], PrRaw)
# })

# Likelihood of Rootnode -------------------------------------------------------
# test_that("Log-Likelihood", {
  ll0 <- LogLike(O, psi = psi, mu_d = mu, mu_s = mu, eta = eta, Pi = Pi)$ll
  
  PI   <- aphylo:::root_node_prob(Pi, S)
  root <- O$pseq[length(O$pseq)]
  ll1  <- log(sum(Pr[root, , drop = TRUE] * PI))
  
  expect_equal(abs(ll1 - ll0), 0, tol = errtol)
# })

# Proportionally when fixing etas or psi ---------------------------------------
# The baseline model, compared against the one with eta, when fixing
# eta0 = eta1 = 1/2, the likelihoods should be proportional to a constant equal
# to (1/2)^(P*#leafs).

# test_that("Leaf Probabilities", {
  
  # Setting eta to be 1/2
  eta    <- c(.5, .5)
  
  # Checking cases compared to the states (0,0) (1,0) (0,1) (1,1)
  PrRaw   <- Pr
  PrRaw[] <- 1
  
  # Node 3 (0,0)
  PrRaw[4, 1] <- (1 - psi[1])^2 
  PrRaw[4, 2] <- psi[2] * (1 - psi[1]) 
  PrRaw[4, 3] <- (1 - psi[1]) * psi[2] 
  PrRaw[4, 4] <- psi[2]^2 
  
  # Node 4 (0,1)
  PrRaw[5, 1] <- (1 - psi[1])*psi[1] 
  PrRaw[5, 2] <- psi[2] * psi[1] 
  PrRaw[5, 3] <- (1 - psi[1]) * (1 - psi[2]) 
  PrRaw[5, 4] <- psi[2] * (1 - psi[2])
  
  # Node 5 (1,0)
  PrRaw[6, 1] <- psi[1]*(1 - psi[1]) 
  PrRaw[6, 2] <- (1 - psi[2]) * (1 - psi[1]) 
  PrRaw[6, 3] <- psi[1] * psi[2] 
  PrRaw[6, 4] <- (1 - psi[2]) * psi[2] 
  
  # Node 6 (1,1)
  PrRaw[7, 1] <- psi[1] ^ 2 
  PrRaw[7, 2] <- (1 - psi[2]) * psi[1] 
  PrRaw[7, 3] <- psi[1] * (1 - psi[2]) 
  PrRaw[7, 4] <- (1 - psi[2]) ^ 2 
  
  M  <- aphylo:::prob_mat(mu)
  
  # Node 1 (0,0)
  of1 <- 
    PrRaw[4:5,1]*M[1,1]*M[1,1] + # (0,0), (0,0)
    PrRaw[4:5,2]*M[1,2]*M[1,1] + # (0,0), (1,0)
    PrRaw[4:5,3]*M[1,1]*M[1,2] + # (0,0), (0,1)
    PrRaw[4:5,4]*M[1,2]*M[1,2]   # (0,0), (1,1)
  
  PrRaw[2,1] <- prod(of1)
  
  # Node 2 (0,0)
  of1 <- 
    PrRaw[6:7,1]*M[1,1]*M[1,1] + # (0,0), (0,0)
    PrRaw[6:7,2]*M[1,2]*M[1,1] + # (0,0), (1,0)
    PrRaw[6:7,3]*M[1,1]*M[1,2] + # (0,0), (0,1)
    PrRaw[6:7,4]*M[1,2]*M[1,2]   # (0,0), (1,1)
  
  PrRaw[3,1] <- prod(of1)
  
  # Node 1 (1,0)
  of1 <- 
    PrRaw[4:5,1]*M[2,1]*M[1,1] + # (1,0), (0,0)
    PrRaw[4:5,2]*M[2,2]*M[1,1] + # (1,0), (1,0)
    PrRaw[4:5,3]*M[2,1]*M[1,2] + # (1,0), (0,1)
    PrRaw[4:5,4]*M[2,2]*M[1,2]   # (1,0), (1,1)
  
  PrRaw[2,2] <- prod(of1)
  
  # Node 2 (1,0)
  of1 <- 
    PrRaw[6:7,1]*M[2,1]*M[1,1] + # (1,0), (0,0)
    PrRaw[6:7,2]*M[2,2]*M[1,1] + # (1,0), (1,0)
    PrRaw[6:7,3]*M[2,1]*M[1,2] + # (1,0), (0,1)
    PrRaw[6:7,4]*M[2,2]*M[1,2]   # (1,0), (1,1)
  
  PrRaw[3,2] <- prod(of1)
  
  # Node 1 (0,1)
  of1 <- 
    PrRaw[4:5,1]*M[1,1]*M[2,1] + # (0,1), (0,0)
    PrRaw[4:5,2]*M[1,2]*M[2,1] + # (0,1), (1,0)
    PrRaw[4:5,3]*M[1,1]*M[2,2] + # (0,1), (0,1)
    PrRaw[4:5,4]*M[1,2]*M[2,2]   # (0,1), (1,1)
  
  PrRaw[2,3] <- prod(of1)
  
  # Node 2 (0,1)
  of1 <- 
    PrRaw[6:7,1]*M[1,1]*M[2,1] + # (0,1), (0,0)
    PrRaw[6:7,2]*M[1,2]*M[2,1] + # (0,1), (1,0)
    PrRaw[6:7,3]*M[1,1]*M[2,2] + # (0,1), (0,1)
    PrRaw[6:7,4]*M[1,2]*M[2,2]   # (0,1), (1,1)
  
  PrRaw[3,3] <- prod(of1)
  
  # Node 1 (1,1)
  of1 <- 
    PrRaw[4:5,1]*M[2,1]*M[2,1] + # (1,1), (0,0)
    PrRaw[4:5,2]*M[2,2]*M[2,1] + # (1,1), (1,0)
    PrRaw[4:5,3]*M[2,1]*M[2,2] + # (1,1), (0,1)
    PrRaw[4:5,4]*M[2,2]*M[2,2]   # (1,1), (1,1)
  
  PrRaw[2,4] <- prod(of1)
  
  # Node 2 (1,1)
  of1 <- 
    PrRaw[6:7,1]*M[2,1]*M[2,1] + # (1,1), (0,0)
    PrRaw[6:7,2]*M[2,2]*M[2,1] + # (1,1), (1,0)
    PrRaw[6:7,3]*M[2,1]*M[2,2] + # (1,1), (0,1)
    PrRaw[6:7,4]*M[2,2]*M[2,2]   # (1,1), (1,1)
  
  PrRaw[3,4] <- prod(of1)
  
  # Node 0 (0,0)
  of1 <- 
    PrRaw[2:3,1]*M[1,1]*M[1,1] + # (0,0), (0,0)
    PrRaw[2:3,2]*M[1,2]*M[1,1] + # (0,0), (1,0)
    PrRaw[2:3,3]*M[1,1]*M[1,2] + # (0,0), (0,1)
    PrRaw[2:3,4]*M[1,2]*M[1,2]   # (0,0), (1,1)
  
  PrRaw[1,1] <- prod(of1)
  
  # Node 0 (1,0)
  of1 <- 
    PrRaw[2:3,1]*M[2,1]*M[1,1] + # (1,0), (0,0)
    PrRaw[2:3,2]*M[2,2]*M[1,1] + # (1,0), (1,0)
    PrRaw[2:3,3]*M[2,1]*M[1,2] + # (1,0), (0,1)
    PrRaw[2:3,4]*M[2,2]*M[1,2]   # (1,0), (1,1)
  
  PrRaw[1,2] <- prod(of1)
  
  # Node 0 (0,1)
  of1 <- 
    PrRaw[2:3,1]*M[1,1]*M[2,1] + # (0,1), (0,0)
    PrRaw[2:3,2]*M[1,2]*M[2,1] + # (0,1), (1,0)
    PrRaw[2:3,3]*M[1,1]*M[2,2] + # (0,1), (0,1)
    PrRaw[2:3,4]*M[1,2]*M[2,2]   # (0,1), (1,1)
  
  PrRaw[1,3] <- prod(of1)
  
  # Node 0 (1,1)
  of1 <- 
    PrRaw[2:3,1]*M[2,1]*M[2,1] + # (1,1), (0,0)
    PrRaw[2:3,2]*M[2,2]*M[2,1] + # (1,1), (1,0)
    PrRaw[2:3,3]*M[2,1]*M[2,2] + # (1,1), (0,1)
    PrRaw[2:3,4]*M[2,2]*M[2,2]   # (1,1), (1,1)
  
  PrRaw[1,4] <- prod(of1)
  
  # Computing likelihood
  ll0 <- LogLike(O, psi = psi, mu_d = mu, mu_s = mu, eta = eta, Pi = Pi)$ll
  
  PI   <- aphylo:::root_node_prob(Pi, S)
  root <- O$pseq[length(O$pseq)]
  ll1  <- log(sum(PrRaw[1, , drop = TRUE] * PI))
  
  # This likelihoods are proportional to (1/2)^(P*#leafs)
  expect_equal((exp(ll0)/exp(ll1))^(1/(2*4)), 0.5, tol = errtol)
  
# })

# ------------------------------------------------------------------------------
# test_that("Joint loglike for multiAphylo", {
  
set.seed(113)
trees <- lapply(sample(20:50, 10, TRUE), raphylo)
trees <- do.call(c, trees)

psi  <- c(0.01, 0.05)
mu   <- c(0.03, 0.10)
Pi   <- .2
eta  <- c(.7, .9)
ans0 <- LogLike(trees, psi = psi, mu_d = mu, mu_s = mu, Pi = Pi, eta = eta, verb_ans = FALSE)
ans1 <- 0
for (i in seq_along(trees)) {
  
  ans1 <- ans1 +
    LogLike(
      tree     = trees[[i]],
      psi      = psi,
      mu_d     = mu,
      mu_s     = mu,
      Pi       = Pi,
      eta      = eta,
      verb_ans = FALSE
      )$ll
  
}

expect_equal(ans1, ans0$ll)
  
# })

  
# Test for multiphylo_pruner ---------------------------------------------------

trees_pruner <- new_aphylo_pruner(trees)

ans2 <- LogLike(
  trees_pruner, psi = psi, mu_d = mu, mu_s = mu, Pi = Pi, eta = eta,
  verb_ans = FALSE)

expect_equal(ans0$ll, ans2$ll)

