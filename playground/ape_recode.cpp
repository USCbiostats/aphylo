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
  int m = ans.n_rows;
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
x10 <- ape::rtree(10)
x100 <- ape::rtree(100)
x1000 <- ape::rtree(1000)

microbenchmark::microbenchmark(
  `  10` = ape_recode(x10$edge, max(x10$edge)),
  ` 100` = ape_recode(x100$edge, max(x100$edge)),
  `1000` = ape_recode(x1000$edge, max(x1000$edge)),
  unit = "s", times=1e3
) 

*/