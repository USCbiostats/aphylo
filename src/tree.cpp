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

// [[Rcpp::export(name = ".recode_as_po")]]
IntegerMatrix recode_as_po(
    const arma::imat & edges
  ) {

  // Creating saving storage
  unsigned int N = edges.n_rows;
  
  // Original positions
  arma::icolvec positions(edges.n_rows, arma::fill::ones);
  positions = arma::cumsum(positions) - 1;
  
  // Temporary edgelist
  arma::imat edges0 = edges;
  
  IntegerMatrix edges1(N, 2u);
  
  // Fetching labels and computing the indegree vector
  // L    : Labels,
  // Lans : Resulting labels
  arma::ivec L = arma::unique(arma::vectorise(edges0));
  IntegerVector Lans(L.size());  
  arma::uvec nparents = fast_table_using_labels(edges0.col(1u), L);
  edges0 = arma::join_rows(edges0, positions);

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
  else if (counter > 1u)
    Rcpp::stop("There is more than 1 root node. Multiply rooted trees are not supported.");
  
  unsigned int i, j, nleafs = 0u;
  for (i = 0u; i < Lans.size(); i++) {
    
    // Find offsprings: Vector of individuals equal to the label of Lans.at(i)
    arma::uvec offspring = arma::find(edges0.col(0u) == Lans.at(i));
    
    // If no offspring (then leaf), then continue
    if (offspring.n_rows == 0u) {
      nleafs++;
      continue;
    }

    // Add them to the list
    for (j = 0; j < offspring.size(); j++) {

      // Adding the index
      Lans.at(counter++) = edges0(offspring.at(j), 1u);
      
      // Adding the
      edges1.at(edges0.at(offspring.at(j), 2u), 0u) = i;
      edges1.at(edges0.at(offspring.at(j), 2u), 1u) = counter - 1u;

    }
    
    // Removing from E
    j = offspring.size();
    while (j != 0)
      edges0.shed_row(offspring.at(--j));
      
  }

  // Creating nametags
  StringVector labels(Lans.size());
  for (int i = 0; i< (int) Lans.size(); i++) {
    char lab[10];
    sprintf(&(lab[0]), "%i", Lans.at(i));
    labels[i] = lab;
  }
  
  // Returning
  edges1.attr("labels") = labels;
  
  return edges1;
  
}


// [[Rcpp::export(name = ".list_offspring")]]
List list_offspring(IntegerMatrix E, int n) {
  std::vector< std::vector<int> > ans(n);
  
  for (int i = 0; i < E.nrow(); i++)
    ans.at(E.at(i, 0) - 1).push_back(E.at(i, 1));
  
  List O(n);
  for (int i = 0; i < n; i++)
    O.at(i) = Rcpp::wrap(ans.at(i));
  
  return O;
}
