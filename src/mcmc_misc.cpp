#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector normal_prop(
    const NumericVector & x,
    const NumericVector & ub,
    const NumericVector & lb,
    double scale
) {
  
  int K = x.size();
  
  // Proposal
  NumericVector ans = x + rnorm(K, 0, 1)*scale;
  
  // Reflexion adjustment
  for (int k=0; k<K; k++) {
    
    while( (ans[k] > ub[k]) | (ans[k] < lb[k]) ) {
      
      if (ans[k] > ub[k]) {
        ans[k] = 2.0*ub[k] - ans[k];
      } else {
        ans[k] = 2.0*lb[k] - ans[k];
      }  
      
    }
      
  }
    
  
  return ans;
}
