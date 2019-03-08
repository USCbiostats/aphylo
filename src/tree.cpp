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

typedef std::vector< std::vector<int> > stdintvec;

// [[Rcpp::export(name = ".list_offspring")]]
ListOf<IntegerVector> list_offspring(IntegerMatrix E, int n) {
  stdintvec ans(n);
  
  for (int i = 0; i < E.nrow(); i++)
    ans.at(E.at(i, 0) - 1).push_back(E.at(i, 1));
  
  List O(n);
  for (int i = 0; i < n; i++)
    O.at(i) = Rcpp::wrap(ans.at(i));
  
  return wrap(O);
}
