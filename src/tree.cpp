#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
IntegerMatrix fast_table(
    const arma::ivec & x
  ) {

  arma::ivec ids = unique(x);
  IntegerMatrix ans(ids.size(), 2u);
  arma::imat x0(x.size(), 1u);
  x0.col(0u) = x;
  
  unsigned int N = ids.size(), i, j;
  
  for (i = 0u; i < N; i++) {
    
    // Filling the first row of the output
    ans.at(i,0u) = ids.at(i);
    ans.at(i,1u) = 0;
    
    // Looping through xsize
    j = 0u;
    while (j < x0.n_rows) {
      
      // If id of i is in x0, then add it!
      if (x0.at(j, 0u) == ids.at(i)) {
        
        // Incrementing counter and removing the row
        ++ans.at(i, 1u);
        x0.shed_row(j);
        
      } else j++;
    }
  }
  
  return ans;
}

// [[Rcpp::export]]
arma::uvec fast_table_using_labels(
    const arma::ivec & x,
    const arma::ivec & ids
) {
  
  arma::imat x0(x.size(), 1u);
  arma::uvec ans(ids.size());
  x0.col(0u) = x;
  
  unsigned int N = ids.size(), i, j;
  
  for (i = 0u; i < N; i++) {
    
    // Filling the first row of the output
    ans.at(i) = 0;
    
    // Looping through xsize
    j = 0u;
    while (j < x0.n_rows) {
      
      // If id of i is in x0, then add it!
      if (x0.at(j, 0u) == ids.at(i)) {
        
        // Incrementing counter and removing the row
        // that we just counted
        ++ans.at(i);
        x0.shed_row(j);
        
      } else ++j;
    }
  }
  
  return ans;
}

//' Recodes an edgelist as a Partially Ordered Tree
//' 
//' The function \code{\link{new_aphylo}} uses this function to make sure that
//' the edgelist that is passed makes a partial order. This is a requirement
//' for the peeling algorithm, which is used explicitly in the \code{LogLike}
//' function.
//' 
//' @details
//' The recoded edgelist is such that in all rows the first element, parent
//' node, as a label that is less than the second element, the offspring, a
//' partial order.
//' 
//' @return A matrix of the same dimension as \code{edges}, an edgelist, recoded
//' to form a partial order. Besides of been of class \code{matrix}, the resulting
//' object is also of class \code{po_tree} and has an aditional attribute:
//' \item{labels}{Named integer vector of size n. Original labels of the edgelist
//' where the names are from 0 to \code{n}.}
//' 
//' @template parameters
//' @templateVar edges 1
//' @export
//' @family Data management functions
//' @examples
//' # Recoding an ape tree -----------------------------------------------------
//' 
//' set.seed(1122233)
//' apetree <- ape::rtree(5)
//' potree  <- as_po_tree(apetree$edge)
//' 
//' apetree$edge
//' potree
//' 
//' # Going back
//' potree[] <- attr(potree, "labels")[potree[] + 1]
//' potree # Ordering is a little off, but is the same tree
//' 
// [[Rcpp::export]]
IntegerMatrix as_po_tree(
    const arma::imat & edges
  ) {

  
  // Creating saving storage
  unsigned int N = edges.n_rows;
  arma::imat edges0(edges);
  
  IntegerMatrix edges1(N, 2u);
  
  // Fetching labels and computing the indegree vector
  arma::ivec L = arma::unique(arma::vectorise(edges0));
  IntegerVector Lans(L.size());  
  arma::uvec nparents = fast_table_using_labels(edges0.col(1u), L);
  

  // Tagging possible root nodes and computing number of parents
  unsigned int counter = 0u;
  for (unsigned int i = 0u; i < L.size(); i++)
    if      (nparents.at(i) == 0u) Lans.at(counter++) = L.at(i);
    else if (nparents.at(i) == 1u) continue;
    else stop("Node with label: %i has %i parents. Multiple parents is not supported yet.",
              L.at(i), nparents.at(i));
  
  // Checking if all this makes sense
  if (counter == 0u)
    Rcpp::stop("There are no root nodes (i.e. indegree == 0).");
  
  unsigned int i, j, iE = 0u;
  for (i = 0u; i < Lans.size(); i++) {
    
    // Find offsprings: Vector of individuals equal to 
    arma::uvec offspring = arma::find(edges0.col(0u) == Lans.at(i));
    
    // If no offspring, then continue
    if (offspring.n_rows == 0u)
      continue;
    
    // Add them to the list
    j = offspring.size();
    while (j != 0) {

      // Adding the index
      Lans.at(counter++) = edges0(offspring.at(--j), 1u);
      
      // Adding the
      edges1.at(iE,   0u)   = i;
      edges1.at(iE++, 1u) = counter - 1u;

      // Removing it from E
      edges0.shed_row(offspring.at(j));

    }
      
  }
  
  // Creating nametags
  StringVector names(Lans.size());
  for (int i = 0; i< (int) Lans.size(); i++) {
    char name[10];
    sprintf(&(name[0]), "%i", i);
    names[i] = name;
  }
  
  Lans.attr("names") = names;
  
  // Returning
  edges1.attr("labels") = Lans;
  edges1.attr("class")  = CharacterVector::create(
    "po_tree", "matrix"
  );
  
  return edges1;
  
}


// [[Rcpp::export]]
List list_offspring(
    const arma::umat & edges
) {
  
  unsigned int n = edges.max() - edges.min() + 1u, i;
  std::vector< std::vector< arma::uword > > offspring(n);
  
  // Listing offsprings
  for (i = 0u; i < edges.n_rows; i++) {
    // Adding the offspring
    offspring.at(edges.at(i, 0u)).push_back(edges.at(i, 1u));
  }
  
  // Coercing into a list
  List ans(n);
  for (i = 0u; i < n; i++) {
    if (offspring.at(i).size() == 0) ans.at(i) = R_NilValue;
    else ans.at(i) = arma::conv_to< arma::urowvec >::from( offspring.at(i) );
  }
  
  return ans;
}


