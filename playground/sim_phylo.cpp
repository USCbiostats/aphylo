// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::umat sim_phylo(
    const List       & offspring,
    const arma::ivec & noffspring,
    const arma::vec  & psi,
    const arma::vec  & mu,
    const arma::vec  & Pi
) {
  
  
  // Vessels
  int N = offspring.size(), N_o;
  arma::umat ans(N,1);
  ans.fill(9u);
  
  for (int p=0; p<1; p++) {
    // Root node function
    ans.at(0,p) = (Pi.at(0) > unif_rand())? 1u : 0u;
    
    // Assigning probabilities to their offspring
    for (int i=0; i<N; i++) {
      
      // Leaf nodes have no offspring. So this is when we include the miss
      // classification factor
      N_o = noffspring.at(i);
      if (!N_o) {
        
        // Gain Probability
        if ((ans.at(i, p) == 0u) && (psi.at(0) > unif_rand()) ) 
          ans.at(i, p) = 1u;
        // Loss probability
        else if ((ans.at(i, p) == 1u) && (psi.at(1) > unif_rand()) ) 
          ans.at(i, p) = 0u;
          
        continue;
      }
        
      // Getting the list of offspring
      arma::ivec O = offspring.at(i);
      
      // Looping through offspring
      for (int o=0; o<N_o; o++) {
        
        // If there is a function
        if (ans.at(i,p) == 1u) 
          // Loss probabilities
          ans.at(O.at(o),p) = (mu.at(1) > unif_rand())? 0u : 1u;
        else 
          // Gain Probabilities
          ans.at(O.at(o),p) = (mu.at(0) > unif_rand())? 1u : 0u;
      }
    }
  }
  
  return ans;
}

/***R
library(phylogenetic)

data(experiment)
data(tree)

# Preprocessing the data
O <- get_offspring(
  experiment, "LeafId", 
  tree, "NodeId", "ParentId"
)

ans <- with(O, sim_phylo(offspring, noffspring, psi = c(.001,.05)*0, mu = c(.5, .5), Pi = c(1,0)));table(ans)

*/