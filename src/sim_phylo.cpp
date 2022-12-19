#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(name=".sim_fun_on_tree")]]
IntegerMatrix sim_fun_on_tree(
    const List       & offspring,
    const IntegerVector & types, 
    const IntegerVector & pseq,
    const NumericVector  & psi,
    const NumericVector  & mu_d,
    const NumericVector  & mu_s,
    const NumericVector  & eta,
    const NumericVector  & Pi,
    int P = 1
) {
  
  if (P > 9999)
    stop("This is nuts, simulating a tree with more than 9999 cannot be done.");
  
  // Vessels
  int N = offspring.size(), N_o;
  std::vector< const NumericVector* > mu(2u);
  mu[0] = & mu_d;
  mu[1] = & mu_s;
  IntegerMatrix ans(N,P);
  ans.fill(9u);
  
  typedef IntegerVector::const_iterator iviter;
  for (int p=0; p<P; p++) {
    
    // Root node function
    ans.at(pseq.at(0) - 1, p) = (Pi.at(0) > unif_rand())? 1u : 0u;
    
    // Assigning probabilities to their offspring.
    for (iviter i = pseq.begin(); i != pseq.end(); ++i) {
      
      // Leaf nodes have no offspring. So this is when we include the miss
      // classification factor
      N_o = Rf_length(offspring.at(*i - 1u));
      if (!N_o) {
        
        if (ans.at(*i - 1u, p) == 0u) {
          
          // Mislabelling a zero
          ans.at(*i - 1u, p) = (psi.at(0) > unif_rand()) ? 1u : 0u;
          
        } else if (ans.at(*i - 1, p) == 1u)  {
          
          // Mislabelling a one
          ans.at(*i - 1u, p) = (psi.at(1u) > unif_rand()) ? 0u : 1u;
          
        } else
          stop("Skipping a leaf node.");
          
        // Likelihood of not annotating
        if (eta.at(ans.at(*i - 1u, p)) < unif_rand()) 
          ans.at(*i - 1u, p) = 9u;
        
        continue;
      }
        
      // Getting the list of offspring
      IntegerVector O = offspring.at(*i - 1);
      
      // Looping through offspring
      for (iviter o = O.begin(); o != O.end(); o++) {
        
        // If there is a function
        if (ans.at(*i - 1u, p) == 1u)      // Loss probabilities
          ans.at(*o - 1u, p) = (mu[types[*i - 1u]]->at(1u) > unif_rand())? 0u : 1u; 
        else if (ans.at(*i - 1, p) == 0u) // Gain Probabilities
          ans.at(*o - 1u, p) = (mu[types[*i - 1u]]->at(0u) > unif_rand())? 1u : 0u; 
        else
          stop("Skipping an internal node.");
      }
    }
  }
  
  // Creating a nice set of names
  StringVector fnames(P, "fun");
  for (int p = 0; p<P ; p++) {
    
    size_t pow10 = 4 - ((p <= 10) ? 1 : static_cast<int>(
      std::ceil(std::log10(p))
      ));

    // Pasting zeros
    for (size_t pp = 0; pp < pow10; ++pp)
      fnames[p] += std::string("0");
    
    fnames[p] += std::to_string(p);
    
  }
  
  ans.attr("dimnames") = List::create(
    R_NilValue,
    fnames
    );
  
  // Needed to clean up afterwards
  mu[0u] = nullptr;
  mu[1u] = nullptr;
  
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
    if (j != (int) (N.size() - 1u))
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

