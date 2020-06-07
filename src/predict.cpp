#include <Rcpp.h>
#include "misc.h"
using namespace Rcpp;

// [[Rcpp::export(name = ".posterior_prob", rng=false)]]
List posterior_prob(
    const NumericMatrix  & Pr_postorder,
    const std::vector< unsigned int > & types,
    const NumericVector  & mu_d,
    const NumericVector  & mu_s,
    const double     & Pi,
    const IntegerVector & pseq,
    const List       & offspring
) {
  
  // Creating output objects
  NumericMatrix Pr_preorder(Pr_postorder.nrow(), Pr_postorder.ncol());
  for (auto iter = Pr_preorder.begin(); iter != Pr_preorder.end(); ++iter)
    *iter = 0.0;
  
  IntegerVector preorder(pseq.size());
  NumericVector Posterior(preorder.size());
  
  // Creating Matrix of probabilities
  std::vector< const NumericVector * > mu(2);
  mu[0] = & mu_d;
  mu[1] = & mu_s;
  
  // Generating the preorder
  std::reverse_copy(pseq.begin(), pseq.end(), preorder.begin());
  
  // Iterators definitions
  typedef IntegerVector::const_iterator iviter;
  typedef IntegerVector::const_iterator Riviter;
  
  // Starting with the root node
  Pr_preorder.at(preorder.at(0u) - 1u, 0u) = 
    Pr_postorder.at(preorder.at(0u) - 1u, 0u)*(1.0 - Pi);
  
  Pr_preorder.at(preorder.at(0u) - 1u, 1u) = 
    Pr_postorder.at(preorder.at(0u) - 1u, 1u)*Pi;
  
  // Computing posterior probabilities
  Posterior.at(preorder.at(0u) - 1u) = Pr_preorder.at(preorder.at(0u) - 1u, 1u)/(
    Pr_preorder.at(preorder.at(0u) - 1u, 1u) +
    Pr_preorder.at(preorder.at(0u) - 1u, 0u)
  );
  
  double D_n_complement_x_n;
  
  for (iviter n = preorder.begin(); n != preorder.end(); ++n) {
    
    // If no offsprings, then nothing to do
    if (! (bool) Rf_length(offspring.at(*n - 1u))) {
      continue;
    }
    
    // Obtaining list of offspring <- this can be improved (speed)
    // can create an std vector of size n
    IntegerVector O(offspring.at(*n - 1u));
    
    for (Riviter o = O.begin(); o != O.end(); ++o) {
      
      // Computing the joint (D_n^c, x_n = 0)
      D_n_complement_x_n = 
        Pr_preorder.at(*n - 1u, 0) / 
        (
            Pr_postorder.at(*o - 1u, 1)*mu[types[*n - 1u]]->at(0) + 
              Pr_postorder.at(*o - 1u, 0)*(1.0 - mu[types[*n - 1u]]->at(0))
        ) * (1 -mu[types[*n - 1u]]->at(0)) +
          Pr_preorder.at(*n - 1u, 1) / 
          (
              Pr_postorder.at(*o - 1u, 1)*(1.0-mu[types[*n - 1u]]->at(1)) + 
                Pr_postorder.at(*o - 1u, 0)*mu[types[*n - 1u]]->at(1)
          ) * mu[types[*n - 1u]]->at(1)
      ;
      
      
      Pr_preorder.at(*o - 1u, 0) = Pr_postorder.at(*o - 1u, 0)*
        D_n_complement_x_n;
      
      // Computing the joint (D_n^c, x_n = 1)
      D_n_complement_x_n = 
        Pr_preorder.at(*n - 1u, 0) / 
        (
            Pr_postorder.at(*o - 1u, 1)*mu[types[*n - 1u]]->at(0) + 
              Pr_postorder.at(*o - 1u, 0)*(1.0 - mu[types[*n - 1u]]->at(0))
        ) * mu[types[*n - 1u]]->at(0) +
          Pr_preorder.at(*n - 1u, 1) / 
          (
              Pr_postorder.at(*o - 1u, 1)*(1.0-mu[types[*n - 1u]]->at(1)) + 
                Pr_postorder.at(*o - 1u, 0)*mu[types[*n - 1u]]->at(1)
          ) * (1 - mu[types[*n - 1u]]->at(1))
        ;
      
      
      Pr_preorder.at(*o - 1u, 1) = Pr_postorder.at(*o - 1u, 1)*
        D_n_complement_x_n;
      
      
      // Computing posterior probabilities
      Posterior.at(*o - 1) = Pr_preorder.at(*o - 1u, 1)/
        (Pr_preorder.at(*o - 1u, 1) + Pr_preorder.at(*o - 1u, 0));
      
    }
    
  }
  
  
  return List::create(
    _["joint_prob"] = Pr_preorder,
    _["posterior"]  = Posterior,
    _["pseq"]       = preorder
  );
  
}
