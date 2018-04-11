#' Posterior probabilities based on parameter estimates
#' @param atree A tree of class [aphylo]
#' @template parameters
#' @templateVar psi 1
#' @templateVar mu 1
#' @templateVar Pi 1
#' @templateVar eta 1
#' @name predict-bis
NULL

#' @rdname predict-bis
#' @export
predict_pre_order <- function(atree, psi, mu, eta, Pi) {
  
  l <- LogLike(atree, psi, mu, eta, Pi)
  
  # const arma::mat  & Pr_postorder,
  # const arma::vec  & mu,
  # const double     & Pi,
  # const arma::ivec & pseq,
  # const List       & offspring
  
  .posterior_prob(l$Pr, mu, Pi, atree$pseq, atree$offspring)
  
}

#' @rdname predict-bis
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
  PSI <- aphylo:::prob_mat(psi)
  MU  <- aphylo:::prob_mat(mu)
  
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