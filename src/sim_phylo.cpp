#include <RcppArmadillo.h>
using namespace Rcpp;

int sample_int(int n) {
  return (int) floor(unif_rand() * (double) n);
}

//' Simulate functions on a ginven tree
//' 
//' @param offspring A List of length \eqn{N} with the set of offspring of
//' each node.
//' @param psi A numeric vector of length 2 (see details).
//' @param mu A numeric vector of length 2 (see details).
//' @param Pi A numeric vector of length 2 (see details).
//' @param P Integer scalar. Number of functions to simulate.
//' 
//' @details
//' 
//' Using the model described in the vignette
//' \href{../doc/peeling_phylo.html}{peeling_phylo.html}
//' 
//' @return An matrix of size \code{length(offspring)*P} with values 9, 0 and 1
//' indicating \code{"no information"}, \code{"no function"} and \code{"function"}.
//' 
//' @examples
//' # Example 1 ----------------------------------------------------------------
//' # We need to simulate a tree
//' set.seed(1231)
//' newtree <- sim_tree(1e3)
//' 
//' # Preprocessing the data
//'     
//' # Simulating
//' ans <- sim_fun_on_tree(
//'   attr(newtree, "offspring"),
//'   psi = c(.001, .05),
//'   mu = c(.01, .05),
//'   Pi = c(.5, .5)
//'   )
//'       
//' # Tabulating results
//' table(ans)
//' 
//' 
//' @name sim_fun_on_tree
// [[Rcpp::export(name=".sim_fun_on_tree")]]
IntegerMatrix sim_fun_on_tree(
    const List       & offspring,
    const arma::ivec & pseq,
    const arma::vec  & psi,
    const arma::vec  & mu,
    const arma::vec  & Pi,
    int P = 1
) {
  
  
  // Vessels
  int N = offspring.size(), N_o;
  IntegerMatrix ans(N,P);
  ans.fill(9u);
  
  typedef arma::ivec::const_iterator iviter;
  for (int p=0; p<P; p++) {
    
    // Root node function
    ans.at(pseq.at(0) - 1u, p) = (Pi.at(0) > unif_rand())? 1u : 0u;
    
    // Assigning probabilities to their offspring
    for (iviter i = pseq.begin(); i != pseq.end(); i++) {
      
      // Leaf nodes have no offspring. So this is when we include the miss
      // classification factor
      N_o = Rf_length(offspring.at(*i - 1));
      if (!N_o) {
        
        // Gain Probability
        if ((ans.at(*i  - 1, p) == 0u) && (psi.at(0) > unif_rand()) ) 
          ans.at(*i - 1, p) = 1u;
        // Loss probability
        else if ((ans.at(*i - 1, p) == 1u) && (psi.at(1) > unif_rand()) ) 
          ans.at(*i - 1, p) = 0u;
          
        continue;
      }
        
      // Getting the list of offspring
      arma::ivec O = offspring.at(*i - 1);
      
      // Looping through offspring
      for (iviter o = O.begin(); o != O.end(); o++) {
        
        // If there is a function
        if (ans.at(*i - 1, p) == 1u) 
          // Loss probabilities
          ans.at(*o - 1, p) = (mu.at(1) > unif_rand())? 0u : 1u;
        else 
          // Gain Probabilities
          ans.at(*o - 1, p) = (mu.at(0) > unif_rand())? 1u : 0u;
      }
    }
  }
  
  // Creating a nice set of names
  StringVector fnames(P);
  for (int p = 0; p<P ; p++) {
    char name[10];
    sprintf(&(name[0]), "fun%04i", p);
    fnames[p] = name;
  }
  
  ans.attr("dimnames") = List::create(
    R_NilValue,
    fnames
    );
  
  return ans;
}

/***R
library(aphylo)

# Preprocessing the data
tree <- sim_tree(1000)
  
# Simulating
ans <- sim_fun_on_tree(
  attr(tree,"offspring"),
  psi = c(.001, .05),
  mu = c(.05, .005),
  Pi = c(.5, .5)
  )

# Tabulating results
table(ans)

*/



// [[Rcpp::export(name = ".sim_tree")]]
IntegerMatrix sim_tree(int n) {
  
  // Initializing
  std::vector< int > offspring, parent;
  std::vector< int > left(n);
  
  std::vector< std::vector<int> > offspring_list(n*2 - 1);
  
  
  // Filling the temporary ids
  for (int i=0;i<n;i++)
    left.at(i) = n - 1 + i;
  
  int i, j;
  int m = n - 1;
  while (left.size() > 1) {
    
    // Random pick of two
    i = sample_int(left.size());
    j = sample_int(left.size() - 1);
    
    // Decreasing parent node
    --m;
    
    // Adding to the edgelist
    offspring.push_back(left.at(i));
    offspring_list.at(m).push_back(left.at(i));
    left.erase(left.begin() + i);
    
    offspring.push_back(left.at(j));
    offspring_list.at(m).push_back(left.at(j));
    
    parent.push_back(m);
    parent.push_back(m);
    
    left.at(j) = m;
    
  }
  
  // Coercing into a Integer Matrix
  IntegerMatrix edges(offspring.size(),2);
  for (unsigned int i=0; i<offspring.size(); i++) {
    edges.at(i,0) = parent.at(i); // - m;
    edges.at(i,1) = offspring.at(i); // - m;
  }
  
  // Coercing into list
  List O(offspring_list.size());
  for (unsigned int i = 0; i<offspring_list.size(); i++) {
    
    if (!offspring_list.at(i).size())
      O.at(i) = IntegerVector::create();
    else
      O.at(i) = wrap(offspring_list.at(i));
  
  }
  
  // And the vector of labels
  StringVector nnames(O.size());
  for (unsigned int i = 0; i<O.size() ; i++) {
    char name[10];
    sprintf(&(name[0]), "%i", i);
    nnames[i] = name;
  }
  
  edges.attr("labels") = nnames;
  edges.attr("offspring") = O;
  
  return edges;
  
}

/***R

microbenchmark::microbenchmark(
  ape = rtree(1e3),
  phy = sim_tree(1e3),
  unit = "relative", times=1e3
)

*/

