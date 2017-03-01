#include <RcppArmadillo.h>
using namespace Rcpp;

//' Simulate functions on a ginven tree
//' 
//' @param offspring A List of length \eqn{N} with the set of offspring of
//' each node.
//' @param noffspring An integer vector of length \eqn{N} with the number of
//' offspring per node.
//' @param psi A numeric vector of length 2 (see details).
//' @param mu A numeric vector of length 2 (see details).
//' @param Pi A numeric vector of length 2 (see details).
//' @param P Integer scalar. Number of functions to simulate.
//' 
//' @details
//' 
//' Using the model described in the vignette
//' \link[=../doc/peeling_phylo.html]{peeling_phylo.html}
//' 
//' @return An matrix of size \code{length(offspring)*P} with values 9, 0 and 1
//' indicating \code{"no information"}, \code{"no function"} and \code{"function"}.
//' 
//' @export
//' @examples
//' # Example 1 ----------------------------------------------------------------
//' # Loading the data
//' data(experiment)
//' data(tree)
//'   
//' # Preprocessing the data
//' O <- get_offspring(experiment, "LeafId", tree, "NodeId", "ParentId")
//'     
//' # Simulating
//' ans <-
//'   with(O, sim_fun_on_tree(
//'       offspring,
//'       noffspring,
//'       psi = c(.001, .05) * 0,
//'       mu = c(.01, .05),
//'       Pi = c(.5, .5)
//'   ))
//'       
//' # Tabulating results
//' table(ans)
//' 
//' 
// [[Rcpp::export]]
arma::umat sim_fun_on_tree(
    const List       & offspring,
    const arma::ivec & noffspring,
    const arma::vec  & psi,
    const arma::vec  & mu,
    const arma::vec  & Pi,
    int P = 1
) {
  
  
  // Vessels
  int N = offspring.size(), N_o;
  arma::umat ans(N,P);
  ans.fill(9u);
  
  for (int p=0; p<P; p++) {
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
O <- get_offspring(experiment, "LeafId", tree, "NodeId", "ParentId")
  
# Simulating
ans <-
with(O, sim_phylo(
    offspring,
    noffspring,
    psi = c(.001, .05) * 0,
    mu = c(.5, .5),
    Pi = c(1, 0)
))

# Tabulating results
table(ans)

*/