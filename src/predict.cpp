#include <RcppArmadillo.h>
#include "misc.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(name = ".posterior_prob")]]
List posterior_prob(
    const arma::mat  & Pr_postorder,
    const arma::vec  & mu,
    const double     & Pi,
    const arma::ivec & pseq,
    const List       & offspring
) {
  
  // Creating output objects
  arma::mat Pr_preorder(Pr_postorder.n_rows, Pr_postorder.n_cols, arma::fill::zeros);
  arma::ivec preorder(pseq.n_elem);
  arma::vec Posterior(preorder.n_elem);
  
  // Creating Matrix of probabilities
  arma::mat M   = prob_mat(mu);
  
  // Generating the preorder
  std::reverse_copy(pseq.begin(), pseq.end(), preorder.begin());
  
  // Iterators definitions
  typedef arma::ivec::const_iterator iviter;
  typedef IntegerVector::const_iterator Riviter;
  
  // Starting with the root node
  Pr_preorder.at(preorder.at(0u) - 1u, 0u) = 
    Pr_postorder.at(preorder.at(0u) - 1u, 0u)*(1.0 - Pi);
  
  Pr_preorder.at(preorder.at(0u) - 1u, 1u) = 
    Pr_postorder.at(preorder.at(0u) - 1u, 1u)*Pi;
  
  // Computing posterior probabilities
  Posterior.at(preorder.at(0u) - 1) = Pr_preorder.at(preorder.at(0u) - 1u, 1)/(
    Pr_preorder.at(preorder.at(0u) - 1u, 1) +
    Pr_preorder.at(preorder.at(0u) - 1u, 0)
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
            Pr_postorder.at(*o - 1u, 1)*mu.at(0) + 
              Pr_postorder.at(*o - 1u, 0)*(1.0 - mu.at(0))
        ) * (1 -mu.at(0)) +
          Pr_preorder.at(*n - 1u, 1) / 
          (
              Pr_postorder.at(*o - 1u, 1)*(1.0-mu.at(1)) + 
                Pr_postorder.at(*o - 1u, 0)*mu.at(1)
          ) * mu.at(1)
      ;
      
      
      Pr_preorder.at(*o - 1u, 0) = Pr_postorder.at(*o - 1u, 0)*
        D_n_complement_x_n;
      
      // Computing the joint (D_n^c, x_n = 1)
      D_n_complement_x_n = 
        Pr_preorder.at(*n - 1u, 0) / 
        (
            Pr_postorder.at(*o - 1u, 1)*mu.at(0) + 
              Pr_postorder.at(*o - 1u, 0)*(1.0 - mu.at(0))
        ) * mu.at(0) +
          Pr_preorder.at(*n - 1u, 1) / 
          (
              Pr_postorder.at(*o - 1u, 1)*(1.0-mu.at(1)) + 
                Pr_postorder.at(*o - 1u, 0)*mu.at(1)
          ) * (1 - mu.at(1))
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
