// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#ifndef APHYLO_MISC_
#define APHYLO_MISC_

using namespace Rcpp;

arma::umat approx_geodesic(
    const arma::umat & edges,
    unsigned int nsteps = 1e3,
    bool undirected = true,
    bool warn = false
);

arma::mat prob_mat(
    const arma::vec & pr
);

arma::imat states(
    int P
);

#endif
