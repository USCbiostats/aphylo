// g++ -I/usr/local/include  -I"/home/vegayon/phylogenetic/src" -o peeling_phylogenies peeling_phylogenies.cpp

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// What does memset do in the context of line 94
// [[Rcpp::export]]
void mymemset() {
  int N = 5u; 
  int Noff[N];
  
  for (int i=0;i<N;i++)
    Rprintf("Noff[%d]: %d\n", i, Noff[i]);
  
  
  Rprintf("After\n");
  memset(Noff, 0, sizeof(Noff));
  for (int i=0;i<N;i++)
    Rprintf("Noff[%d]: %d\n", i, Noff[i]);
  
  return;
}

/***R
mymemset()
# > mymemset()
# Noff[0]: 150939128
# Noff[1]: 0
# Noff[2]: -2104707307
# Noff[3]: 32592
# Noff[4]: 0
# After
# Noff[0]: 0
# Noff[1]: 0
# Noff[2]: 0
# Noff[3]: 0
# Noff[4]: 0
*/

// [[Rcpp::export]]
List getpeelingseq(IntegerVector id, IntegerVector parent) {
  int N = id.length();
  std::vector< IntegerVector > tmp(N);
  
  // Listing offsprings by parent
  for (int i=0; i<N; i++) {
    if (parent.at(i) == NA_INTEGER) continue;
    tmp.at(parent.at(i)).push_back(i);
  }
  
  // Coercing into a List
  List ans(N);
  for (int i=0; i<N; i++)
    ans.at(i) = tmp.at(i);
  
  return ans;
}


// decomposes index of joint function states xx into vector of single-function
// states x[FF]  
// [[Rcpp::export]]
void GetX(int xx) {
  int F = 3, x[3u];
  for (int f=0; f<F; f++) {
    x[f] = xx%2;
    xx /= 2;
  }	
  
  for (int f=0; f<F; f++)
    Rprintf("x[%d] = %d\n", f, x[f]);
  
  return;
  
}

/***R
GetX(0)
GetX(1)
GetX(2)
GetX(3)
*/