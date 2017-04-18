#include <RcppArmadillo.h>
using namespace Rcpp;

int sample_int(int n) {
  return (int) floor(unif_rand() * (double) n);
}

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
//' \href{../doc/peeling_phylo.html}{peeling_phylo.html}
//' 
//' @return An matrix of size \code{length(offspring)*P} with values 9, 0 and 1
//' indicating \code{"no information"}, \code{"no function"} and \code{"function"}.
//' 
//' @export
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
//'   attr(newtree, "noffspring"),
//'   psi = c(.001, .05),
//'   mu = c(.01, .05),
//'   Pi = c(.5, .5)
//'   )
//'       
//' # Tabulating results
//' table(ans)
//' 
//' 
// [[Rcpp::export]]
IntegerMatrix sim_fun_on_tree(
    const List       & offspring,
    const arma::ivec & noffspring,
    const arma::vec  & psi,
    const arma::vec  & mu,
    const arma::vec  & Pi,
    int P = 1
) {
  
  
  // Vessels
  int N = offspring.size(), N_o;
  IntegerMatrix ans(N,P);
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
  attr(tree,"noffspring"),
  psi = c(.001, .05),
  mu = c(.05, .005),
  Pi = c(.5, .5)
  )

# Tabulating results
table(ans)

*/


//' Random tree generation
//' 
//' By randomly choosing pairs of vertices, this function generates random
//' trees from bottom to top such that the labels of the nodes follow a 
//' partial order in their parent-offspring relation, with the parent
//' always having a lower idlabel than the offspring.
//' 
//' @param n Integer scalar. Number of leaf nodes.
//' 
//' @details The algorithm was implemented as follows
//' 
//' \enumerate{
//'   \item Initialize \code{left =[n*2 - 2,...,(n-1)]} and \code{m = n*2 - 2}, and
//'         initialize the vectors \code{parent} and \code{offspring} to be
//'         empty.
//'   \item While \code{length(left) > 1} do:
//'   \enumerate{
//'     \item Randomly choose a pair \code{(i, j)} from \code{left}
//'     \item Add \code{leaf(i)}, \code{leaf(j)} to the tail of \code{offspring},
//'     \item Decrease \code{m} by 1, and add it two times to the tail of
//'     \code{parent}.
//'     \item Remove \code{(i,j)} from \code{leaf} and add \code{m} to its tail.
//'     \item next
//'   }
//'   
//' }
//' 
//' The \code{\link[ape:rtree]{rtree}} function in the \pkg{ape} package is similar,
//' although the big difference is in the way the labels are stablished. This later
//' point is crucial for both \pkg{phylogenetic} and \pkg{ape} as is a key feature
//' in some (most) of its routines.
//' 
//' @return An matrix of size \code{(n*2 - 2)*2}, an edgelist, with \code{n*2-1} nodes.
//' This, a Directed Acyclic Graph (DAG), as classes \code{matrix} and \code{po_tree}.
//' Also, includes the following attributes:
//' 
//' \item{offspring}{A list of size \code{n*2 - 1} listing node ith's offspring if any.}
//' \item{noffspring}{An integer vector of size \code{n*2 - 1} indicating the number of
//' offspring that each node has.}
//' 
//' @examples
//' # A very simple example ----------------------------------------------------
//' set.seed(1223)
//' newtree <- sim_tree(50)
//' 
//' plot(as.apephylo(newtree))
//' 
//' # This is what you would do in igraph --------------------------------------
//' \dontrun{
//' g   <- ans
//' g[] <- as.character(g)
//' g <- igraph::graph_from_edgelist(g)
//' plot(g, layout = igraph::layout_with_sugiyama(g)[[2]])
//' }
//' 
//' # A performance benchmark with ape::rtree ----------------------------------
//' \dontrun{
//' microbenchmark::microbenchmark(
//' ape = rtree(1e3),
//'   phy = sim_tree(1e3),
//' unit = "relative"
//' )
//' # This is what you would get.
//' Unit: relative
//'   expr     min       lq     mean  median       uq      max neval
//'    ape 14.7598 14.30809 14.30013 16.7217 14.32843 4.754106   100
//'    phy  1.0000  1.00000  1.00000  1.0000  1.00000 1.000000   100
//' }
//' @export
// [[Rcpp::export]]
IntegerMatrix sim_tree(int n) {
  
  // Initializing
  std::vector< int > offspring, parent;
  std::vector< int > left(n);
  
  std::vector< std::vector<int> > offspring_list(n*2 - 1);
  IntegerVector noffspring(offspring_list.size(), 0);
  
  
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
    noffspring.at(m)+=2;
    
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
    if (!offspring_list.at(i).size()) continue;
    
    O.at(i) = wrap(arma::conv_to< arma::urowvec >::from(offspring_list.at(i)));
  }
  
  // Creating the answer
  edges.attr("offspring")  = O;
  edges.attr("noffspring") = noffspring;
  edges.attr("class") = CharacterVector::create(
    "po_tree",
    "matrix"
  );
  
  return edges;
  
}

/***R

microbenchmark::microbenchmark(
  ape = rtree(1e3),
  phy = sim_tree(1e3),
  unit = "relative", times=1e3
)

*/

