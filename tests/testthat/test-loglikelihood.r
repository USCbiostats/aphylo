context("LogLikelihood function")

data("faketree")
data("fakeexperiment")

# Parameters
psi    <- c(.01, .02)
mu     <- c(.05, .1)
Pi     <- 1 - .3
errtol <- 1e-15

O  <- new_aphylo(tip.annotation = fakeexperiment[,-1L], tree = faketree)

S  <- states(2)
Pr <- probabilities(
  with(O, rbind(tip.annotation, node.annotation)),
  O$pseq, mu, psi, S, O$offspring)

# Checking Leaf Probabilities --------------------------------------------------

test_that("Leaf Probabilities", {
  # Checking cases compared to the states (0,0) (1,0) (0,1) (1,1)
  PrRaw   <- Pr
  PrRaw[] <- 1
  
  # Node 3 (0,0)
  PrRaw[4, 1] <- (1 - psi[1]) ^ 2
  PrRaw[4, 2] <- psi[2] * (1 - psi[1])
  PrRaw[4, 3] <- (1 - psi[1]) * psi[2]
  PrRaw[4, 4] <- psi[2] ^ 2
  
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
  
  # These should be identical
  expect_equivalent(PrRaw[4:7,], Pr[c(5:7, 1:4),][4:7,])
})


# Checking Internal Probabilities ----------------------------------------------
test_that("Internal Probabilities", {
  M  <- aphylo:::prob_mat(mu)
  # Checking cases compared to the states (0,0) (1,0) (0,1) (1,1)
  
  PrRaw <- Pr
  PrRaw[4:7,] <- 1
  
  # Node 1 (0,0)
  of1 <- 
    Pr[4:5,1]*M[1,1]*M[1,1] + # (0,0), (0,0)
    Pr[4:5,2]*M[1,2]*M[1,1] + # (0,0), (1,0)
    Pr[4:5,3]*M[1,1]*M[1,2] + # (0,0), (0,1)
    Pr[4:5,4]*M[1,2]*M[1,2]   # (0,0), (1,1)
  
  PrRaw[2,1] <- prod(of1)
  
  # Node 2 (0,0)
  of1 <- 
    Pr[6:7,1]*M[1,1]*M[1,1] + # (0,0), (0,0)
    Pr[6:7,2]*M[1,2]*M[1,1] + # (0,0), (1,0)
    Pr[6:7,3]*M[1,1]*M[1,2] + # (0,0), (0,1)
    Pr[6:7,4]*M[1,2]*M[1,2]   # (0,0), (1,1)
  
  PrRaw[3,1] <- prod(of1)
  
  # Node 1 (1,0)
  of1 <- 
    Pr[4:5,1]*M[2,1]*M[1,1] + # (1,0), (0,0)
    Pr[4:5,2]*M[2,2]*M[1,1] + # (1,0), (1,0)
    Pr[4:5,3]*M[2,1]*M[1,2] + # (1,0), (0,1)
    Pr[4:5,4]*M[2,2]*M[1,2]   # (1,0), (1,1)
  
  PrRaw[2,2] <- prod(of1)
  
  # Node 2 (1,0)
  of1 <- 
    Pr[6:7,1]*M[2,1]*M[1,1] + # (1,0), (0,0)
    Pr[6:7,2]*M[2,2]*M[1,1] + # (1,0), (1,0)
    Pr[6:7,3]*M[2,1]*M[1,2] + # (1,0), (0,1)
    Pr[6:7,4]*M[2,2]*M[1,2]   # (1,0), (1,1)
  
  PrRaw[3,2] <- prod(of1)
  
  # Node 1 (0,1)
  of1 <- 
    Pr[4:5,1]*M[1,1]*M[2,1] + # (0,1), (0,0)
    Pr[4:5,2]*M[1,2]*M[2,1] + # (0,1), (1,0)
    Pr[4:5,3]*M[1,1]*M[2,2] + # (0,1), (0,1)
    Pr[4:5,4]*M[1,2]*M[2,2]   # (0,1), (1,1)
  
  PrRaw[2,3] <- prod(of1)
  
  # Node 2 (0,1)
  of1 <- 
    Pr[6:7,1]*M[1,1]*M[2,1] + # (0,1), (0,0)
    Pr[6:7,2]*M[1,2]*M[2,1] + # (0,1), (1,0)
    Pr[6:7,3]*M[1,1]*M[2,2] + # (0,1), (0,1)
    Pr[6:7,4]*M[1,2]*M[2,2]   # (0,1), (1,1)
  
  PrRaw[3,3] <- prod(of1)
  
  # Node 1 (1,1)
  of1 <- 
    Pr[4:5,1]*M[2,1]*M[2,1] + # (1,1), (0,0)
    Pr[4:5,2]*M[2,2]*M[2,1] + # (1,1), (1,0)
    Pr[4:5,3]*M[2,1]*M[2,2] + # (1,1), (0,1)
    Pr[4:5,4]*M[2,2]*M[2,2]   # (1,1), (1,1)
  
  PrRaw[2,4] <- prod(of1)
  
  # Node 2 (1,1)
  of1 <- 
    Pr[6:7,1]*M[2,1]*M[2,1] + # (1,1), (0,0)
    Pr[6:7,2]*M[2,2]*M[2,1] + # (1,1), (1,0)
    Pr[6:7,3]*M[2,1]*M[2,2] + # (1,1), (0,1)
    Pr[6:7,4]*M[2,2]*M[2,2]   # (1,1), (1,1)
  
  PrRaw[3,4] <- prod(of1)
  
  # Node 0 (0,0)
  of1 <- 
    Pr[2:3,1]*M[1,1]*M[1,1] + # (0,0), (0,0)
    Pr[2:3,2]*M[1,2]*M[1,1] + # (0,0), (1,0)
    Pr[2:3,3]*M[1,1]*M[1,2] + # (0,0), (0,1)
    Pr[2:3,4]*M[1,2]*M[1,2]   # (0,0), (1,1)
  
  PrRaw[1,1] <- prod(of1)
  
  # Node 0 (1,0)
  of1 <- 
    Pr[2:3,1]*M[2,1]*M[1,1] + # (1,0), (0,0)
    Pr[2:3,2]*M[2,2]*M[1,1] + # (1,0), (1,0)
    Pr[2:3,3]*M[2,1]*M[1,2] + # (1,0), (0,1)
    Pr[2:3,4]*M[2,2]*M[1,2]   # (1,0), (1,1)
  
  PrRaw[1,2] <- prod(of1)
  
  # Node 0 (0,1)
  of1 <- 
    Pr[2:3,1]*M[1,1]*M[2,1] + # (0,1), (0,0)
    Pr[2:3,2]*M[1,2]*M[2,1] + # (0,1), (1,0)
    Pr[2:3,3]*M[1,1]*M[2,2] + # (0,1), (0,1)
    Pr[2:3,4]*M[1,2]*M[2,2]   # (0,1), (1,1)
  
  PrRaw[1,3] <- prod(of1)
  
  # Node 0 (1,1)
  of1 <- 
    Pr[2:3,1]*M[2,1]*M[2,1] + # (1,1), (0,0)
    Pr[2:3,2]*M[2,2]*M[2,1] + # (1,1), (1,0)
    Pr[2:3,3]*M[2,1]*M[2,2] + # (1,1), (0,1)
    Pr[2:3,4]*M[2,2]*M[2,2]   # (1,1), (1,1)
  
  PrRaw[1,4] <- prod(of1)
  
  expect_equal(Pr[c(5:7, 1:4),], PrRaw)
})

# Likelihood of Rootnode -------------------------------------------------------

test_that("Log-Likelihood", {
  ll0 <- LogLike(O, psi, mu, Pi)$ll
  
  PI  <- aphylo:::root_node_prob(Pi, S)
  ll1 <- log(sum(Pr[1, , drop = TRUE] * PI))
  
  expect_equal(abs(ll1 - ll0), 0, tol = errtol)
})
