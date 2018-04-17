library(ape)
library(aphylo)

# Creating the tree ------------------------------------------------------------
# tree <- matrix(c(3, 1, 3, 2), ncol=2,byrow=TRUE)
# tree <- as.phylo(tree)
# X <- c(1, 0)
tree <- matrix(c(5, 1, 5, 2, 1, 3, 1, 4), ncol=2, byrow = TRUE)
# tree <- matrix(c(5, 6, 5, 1, 6, 2, 6, 3, 6, 4), ncol=2, byrow = TRUE)
tree <- as.phylo(tree)
X    <- c(0, 0, 0)
atree <- new_aphylo(X, tree)

# Model parameters
Pi  <- .3
psi <- c(.01, .02)
mu  <- c(.2, .1)
eta <- c(1, 1)

#' Posterior probabilities based on parameter estimates
#' @param atree A tree of class [aphylo]
#' @template parameters
#' @templateVar psi 1
#' @templateVar mu 1
#' @templateVar Pi
#' @export
predict_brute_force <- function(atree, psi, mu, Pi) {
  
  # Coercing into the true class
  tree <- as.phylo(atree)
  Pi   <- c(1 - Pi, Pi)
  
  # Generating a states matrix
  states <- do.call(expand.grid, rep(list(c(0,1)), ape::Nnode(tree) + ape::Ntip(tree)))
  colnames(states) <- paste0("node", 1:(ncol(states)))
  
  # Computing matrix of probabilities --------------------------------------------
  # 2^(ntips + nnodes) 
  PSI <- prob_mat(psi)
  MU  <- prob_mat(mu)
  
  # For
  Pr <- Pi[states[, ape::Ntip(tree) + 1] + 1]
  for (i in 1:nrow(tree$edge)) {
    e <- tree$edge[i,]
    Pr <- Pr *
      MU[cbind(states[, e[1]], states[, e[2]]) + 1]
  }
  
  # Computing for each possible annotation of the tips
  # 2^Ntip
  Pr <- matrix(rep(Pr, 3), nrow=length(Pr), ncol = 2^ape::Ntip(tree))
  states_tip <- do.call(expand.grid, rep(list(c(0,1)), ape::Ntip(tree)))
  
  for (i in 1:nrow(states_tip)) {
    
    for (j in 1:ape::Ntip(tree)) {
      Pr[, i] <- Pr[, i] *
        PSI[cbind(states[, j] + 1, states_tip[i, j] + 1)]
    }
    
  }
  
  # Computing posterior probabilities
  posterior <- vector("numeric", ncol(states))
  for (i in 1:ncol(states)) {
    state_col <- which(apply(states_tip, 1, function(x) all(x == atree$tip.annotation)))
    state_rows <- which(states[,i] == 1)
    
    posterior[i] <- sum(Pr[state_rows, state_col])/sum(Pr[, state_col])
  }
  
  
  # Results
  list(
    Pr  = Pr,
    row = states,
    col = states_tip,
    posterior = posterior
  )
  
}

ans <- predict_brute_force(atree, psi, mu, Pi)

sum(ans$Pr) # This should add up to 1

# Comparing with aphylo --------------------------------------------------------
l <- LogLike(atree, psi, mu, eta, Pi)
X <- as.vector(atree$tip.annotation)
colSums(ans$Pr)[which(apply(ans$col, 1, function(x) all(x == X)))] -
  exp(l$ll)

exp(l$ll) # Jmorr tree 0.286622598, c(0, 0, 0)

# Conditional probability ------------------------------------------------------

i <- 5

# Matching state
state_col <- which(apply(ans$col, 1, function(x) all(x == X)))
state_rows <- which(ans$row[,i] == 1)

# Pr(root = 0 | data)
sprintf(
  "%.10f",
  1 - sum(ans$Pr[state_rows, state_col])/sum(ans$Pr[, state_col]) # 0.9907319671
)

sprintf(
  "%.10f",
  l$Pr[4,1]*(1-Pi)/
    (l$Pr[4,2]*Pi + l$Pr[4,1]*(1-Pi))
)

i_isone <- which(ans$row[,i] == 1L)
i_iszero <- which(ans$row[,i] == 0L)
X_given_xi1 <- sum(ans$Pr[i_isone, state_col])/sum(ans$Pr[i_isone,])
X_given_xi0 <- sum(ans$Pr[i_iszero, state_col])/sum(ans$Pr[i_iszero,])
xi1 <- sum(ans$Pr[i_isone,])

sprintf(
  "%.10f",
  1 - X_given_xi1/(X_given_xi1 + X_given_xi0*((1-xi1)/xi1))
)




# Jmorr tree
# 
# Pr(z1=0|X)	0.9907319671
# Pr(z2=0|X)	0.9935955962
# Pr(z3=0|X)	0.9931689347
# Pr(z4=0|X)	0.9939582682
# Pr(z5=0|X)	0.9939582682
