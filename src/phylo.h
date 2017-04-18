// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#ifndef APHYLO_PHYLO_
#define APHYLO_PHYLO_

using namespace Rcpp;

arma::imat states(
    int P
  );

arma::mat prob_mat(
    const arma::vec & pr
  );

arma::mat leaf_prob(
    const arma::imat & Z,
    const arma::imat & S,
    const arma::vec  & psi,
    const arma::ivec & noffspring
  );

arma::vec root_node_prob(
    const arma::vec  & pi,
    const arma::imat & S
  );

arma::mat internal_prob(
    arma::mat          Pr,
    const arma::vec  & mu,
    const arma::imat & S,
    const arma::ivec & noffspring,
    const List       & offspring
  );

List LogLike(
    const arma::imat & Z,
    const List       & offspring,
    const arma::ivec & noffspring,
    const arma::vec  & psi,
    const arma::vec  & mu,
    const arma::vec  & Pi,
    bool verb_ans = false
  );


  
#endif
