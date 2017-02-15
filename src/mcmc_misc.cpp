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

// [[Rcpp::export]]
NumericVector gibbs_sampler(
    const Function & fun,
    const NumericVector & x,
    const NumericVector & lb,
    const NumericVector & ub,
    double scale
  ) {
  
  
  int K = x.size();
  
  // Generating proposals
  NumericVector theta0 = x;
  NumericVector theta1 = theta0;
  NumericVector prop = normal_prop(x, ub, lb, scale);
  
  // Which to keep?
  double val0, val1;
  
  for (int k = 0; k < K; k++) {
      theta1.at(k) = prop.at(k);
    
      val0 = as< double >(fun(theta0));
      val1 = as< double >(fun(theta1));
      
      if (exp(val1 - val0) > unif_rand()) {
        theta0.at(k) = prop.at(k);
      } else {
        theta1.at(k) = theta0.at(k);
      }
  }
  
  return theta1;
}
  
  
  
