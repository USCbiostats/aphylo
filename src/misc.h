#include <Rcpp.h>

#ifndef APHYLO_MISC_H
#define APHYLO_MISC_H

using namespace Rcpp;

NumericMatrix prob_mat(
    const NumericVector & pr
);

IntegerMatrix states(
    int P
);

#endif
