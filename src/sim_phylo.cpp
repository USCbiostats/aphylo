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
    ans.at(pseq.at(0) - 1, p) = (Pi.at(0) > unif_rand())? 1u : 0u;
    
    // Assigning probabilities to their offspring.
    for (iviter i = pseq.begin(); i != pseq.end(); i++) {
      
      // Leaf nodes have no offspring. So this is when we include the miss
      // classification factor
      N_o = Rf_length(offspring.at(*i - 1));
      if (!N_o) {
        
        // Gain Probability
        if (ans.at(*i  - 1, p) == 0u) 
          ans.at(*i - 1, p) = (psi.at(0) > unif_rand()) ? 1u : 0u;
        // Loss probability
        else if (ans.at(*i - 1, p) == 1u) 
          ans.at(*i - 1, p) = (psi.at(1) > unif_rand()) ? 0u : 1u;
        else
          stop("Skipping a leaf node.");
          
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
        else if (ans.at(*i - 1, p) == 0u)
          // Gain Probabilities
          ans.at(*o - 1, p) = (mu.at(0) > unif_rand())? 1u : 0u;
        else
          stop("Skipping an internal node.");
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



// [[Rcpp::export(name=".sim_tree")]]
List sim_tree(int n, Function f, bool branches) {
  
  // Initializing set
  std::vector< int > N(n);
  for (int i=0; i<n; i++) {
    N.at(i) = i + 1u;
  }
  
  // Creating labels
  IntegerVector tiplabel = Rcpp::clone(Rcpp::wrap(N));
  IntegerVector nodelabel(n - 1);
  for (int i = 0; i < (n-1); i++)
    nodelabel.at(i) = n + i + 1;
  
  // Initializing edges and pseq
  NumericMatrix E(2*n - 2, 2);
  NumericVector Pseq(2*n - 1);
  
  // Starting the algorithm
  int i, j, p = 0, e = 0;
  int k = n*2 - 1;
  
  while (N.size() > 1) {
    
    // Picking offspring
    i = floor(unif_rand() * N.size());
    j = floor(unif_rand() * (N.size() - 1));
    
    // Adjusting, so don't pick the same
    if (i <= j)
      j++;
    
    // Adding edges
    E.at(e, 0) = k;
    E.at(e++, 1) = N.at(i);
    E.at(e, 0) = k;
    E.at(e++, 1) = N.at(j);
    
    // Adding to peeling seq
    if (N.at(i) <= n)
      Pseq.at(p++) = N.at(i);
    if (N.at(j) <= n)
      Pseq.at(p++) = N.at(j);
    
    Pseq.at(p++) = k;
    
    // Modifying List of nodes:
    // We swap the values so that we can handle the case in which one of them
    // is the last node in N. If that's the case, then we don't need to swap
    // values with the back, which is poped out, otherwise, we need to do so
    // so the last value is kept in the list.
    if (i > j)
      std::swap(i, j);
    if (j != (N.size() - 1))
      N.at(j) = N.back();
    
    N.at(i) = k--;
    
    N.pop_back();
  }
  
  
  // Returning
  List ape;
  if (branches)
    ape = List::create(
      _["edge"]        = E,
      _["edge.length"] = f(2*n - 2),
      _["tip.label"]   = tiplabel,
      _["Nnode"]       = n - 1,
      _["node.label"]  = nodelabel
    );
  else 
    ape = List::create(
      _["edge"]        = E,
      _["tip.label"]   = tiplabel,
      _["Nnode"]       = n - 1,
      _["node.label"]  = nodelabel
    );
  
  
  ape.attr("class") = "phylo";
  ape.attr("order") = "postorder";
  
  return ape;
  
}

/***R

microbenchmark::microbenchmark(
  ape = rtree(1e3),
  phy = sim_tree(1e3),
  unit = "relative", times=1e3
)

*/

