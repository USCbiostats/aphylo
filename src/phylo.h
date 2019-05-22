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

double LogLikei(
    const arma::imat & annotations,
    const ListOf<IntegerVector> & offspring,
    const arma::ivec & pseq,
    const arma::vec  & psi,
    const arma::vec  & mu,
    const arma::vec  & eta,
    double Pi,
    bool check_dims,
    arma::mat & Pr
  );

double LogLike(
    const std::vector< arma::imat > & annotations,
    const std::vector< ListOf<IntegerVector> > & offspring,
    const std::vector< arma::ivec > & pseq,
    const arma::vec & psi,
    const arma::vec & mu,
    const arma::vec & eta,
    double Pi,
    bool verb_ans = false,
    bool check_dims = true
);

  
#endif
