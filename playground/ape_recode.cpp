#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List ape_recode(
  const arma::imat & x,
  int n
) {
  
  // Defining variables
  arma::imat ans;
  ans.copy_size(x);
  int m          = ans.n_rows;
  arma::ivec I(n, arma::fill::zeros);
  
  int o = 0, p = 0;
  // Recoding parents
  for (int i=0; i<m; i++) {
    
    // Setting parent id
    if (I.at(x.at(i, 0) - 1) == 0)
      I.at(x.at(i, 0) - 1) = --p;
    
    ans.at(i, 0) = I.at(x.at(i, 0) - 1);
    
  }

  // Recoding offspring
  for (int i=0; i<m; i++) {
    // Setting offspring id
    if (I[x.at(i, 1) - 1] == 0)
      I[x.at(i, 1) - 1] = ++o;
    
    ans.at(i, 1) = I.at(x.at(i, 1) - 1);
  }
  
  
  // Offspring are alredy from 1 to n, so we need to recode
  // the parents now. 
  for (int i = 0; i<m; i++) {
    if (I.at(x.at(i, 0) - 1) < 0) 
      I.at(x.at(i, 0) - 1) = o - I.at(x.at(i, 0) - 1);
    ans.at(i, 0) = I.at(x.at(i, 0) - 1);
    
    if (I.at(x.at(i, 1) - 1) < 0) 
      I.at(x.at(i, 1) - 1) = o - I.at(x.at(i, 1) - 1);
    ans.at(i, 1) = I.at(x.at(i, 1) - 1);
      
  }
  
  return List::create(
   _["edge"] = ans,
   _["I"] = I
  );
  
  
  
}

/***R
set.seed(1)
x <- ape::rtree(2000)

microbenchmark::microbenchmark(
  ape_recode(x$edge, max(x$edge)), unit = "s"
) 

*/