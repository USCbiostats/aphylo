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

void _fast_table(
    const arma::ivec & x,
    const arma::ivec & ids,
    arma::uvec & ans
) {
  
  arma::imat x0(x.size(), 1u);
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
        ++ans.at(i);
        x0.shed_row(j);
        
      } else j++;
    }
  }
  
  return;
}

//' Recodes an edgelist as a Partially Ordered Tree
//' @template parameters
//' @templateVar edges 1
//' @export
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
  arma::uvec Lans(L.size());  
  arma::uvec nparents(L.size());
  _fast_table(edges0.col(1u), L, nparents);
  

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


