#include <Rcpp.h>
using namespace Rcpp;

int sample_int(int n) {
  return (int) floor(unif_rand() * (double) n);
}

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
//'   \item Initialize \code{left =[0,...,(n-1)]} and \code{m = 0}, and
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
//'   \item Substract \code{m} to each element of both lists, this way the 
//'   root node will have id 0.
//' }
//' 
//' The \code{\link[ape:rtree]{rtree}} function in the \pkg{ape} package is similar,
//' although the big difference is in the way the labels are stablished. This later
//' point is crucial for both \pkg{phylogenetic} and \pkg{ape} as is a key feature
//' in some (most) of its routines.
//' 
//' @return An matrix of size \code{n*2 - 2} with column names \code{"offspring"} and
//' \code{"parent"} representing an edgelist with \code{n*2-1} nodes. A Directed
//' Acyclic Graph (DAG).
//' 
//' @examples
//' # A very simple example ----------------------------------------------------
//' set.seed(1223)
//' ans <- sim_tree(50);range(ans)
//' 
//' # To visualize it, we create a random experimental dataset
//' O <- get_offspring(
//'   data_exper  = data.frame(
//'     f = rep(1L, nrow(ans)),
//'     id = (nrow(ans)-1):nrow(ans)
//'     ),
//'   data_tree   = data.frame(ans),
//'   leafidvar   = "id", 
//'   nodeidvar   = "offspring",
//'   parentidvar = "parent"
//'   )
//' 
//' # Now we plot it as a phylo object
//' library(ape)
//' plot(as.phylo(O))
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
//' # Unit: relative
//' #   expr      min       lq     mean  median       uq      max neval
//' #    ape 33.31827 31.85769 33.04708 30.8261 36.51258 25.75327   100
//' #    phy  1.00000  1.00000  1.00000  1.0000  1.00000  1.00000   100
//' }
//' @export
// [[Rcpp::export]]
IntegerMatrix sim_tree(int n) {
  
  // Initializing
  std::vector< int > source, target;
  std::vector< int > left(n);
  
  // Filling the temporary ids
  for (int i=0;i<n;i++)
    left.at(i) = i;
  
  int i, j;
  int m = 0;
  while (left.size() > 1) {
    
    // Random pick of two
    i = sample_int(left.size());
    j = sample_int(left.size() - 1);
    
    // Adding to the edgelist
    source.push_back(left.at(i));
    left.erase(left.begin() + i);
    
    source.push_back(left.at(j));
    left.at(j) = --m;
    
    target.push_back(m);
    target.push_back(m);
    
  }
  
  // Coercing into a Integer Matrix
  IntegerMatrix ans(source.size(),2);
  for (unsigned int i=0; i<source.size(); i++) {
    ans.at(i,0) = source.at(i) - m;
    ans.at(i,1) = target.at(i) - m;
  }
  
  // Adding some attributes
  ans.attr("dimnames") = List::create(
    R_NilValue,
    CharacterVector::create("offspring","parent")
  );
  
  return ans;
    
}

/***R

microbenchmark::microbenchmark(
  ape = rtree(1e3),
  phy = sim_tree(1e3),
  unit = "relative"
)

*/