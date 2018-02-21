// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#ifndef APHYLO_PHYLO_
#define APHYLO_PHYLO_

using namespace Rcpp;

arma::vec root_node_prob(
    double Pi,
    const arma::imat & S
);

arma::mat probabilities(
    const arma::imat & annotations,
    const arma::ivec & pseq,
    const arma::vec  & psi,
    const arma::vec  & mu,
    const arma::vec  & eta,
    const arma::imat & S,
    const List       & offspring
);

List LogLike(
    const arma::imat & Z,
    const List       & offspring,
    const arma::imat & pseq,
    const arma::vec  & psi,
    const arma::vec  & mu,
    const arma::vec  & eta,
    const arma::vec  & Pi,
    bool verb_ans = false
  );


  
#endif
