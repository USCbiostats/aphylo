#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(rng=false)]]
IntegerMatrix fast_table(
    const IntegerVector & x
  ) {

  IntegerVector ids = unique(x);
  std::sort(ids.begin(), ids.end());
  IntegerMatrix ans(ids.size(), 2u);
  std::vector< int > x0(x.size());
  for (unsigned int i = 0u; i < x0.size(); ++i)
    x0[i] = x.at(i);
  
  unsigned int N = ids.size(), i, j;
  
  for (i = 0u; i < N; i++) {
    
    // Filling the first row of the output
    ans.at(i,0u) = ids.at(i);
    ans.at(i,1u) = 0;
    
    // Looping through xsize
    j = 0u;
    while (j < x0.size()) {
      
      // If id of i is in x0, then add it!
      if (x0[j] == ids.at(i)) {
        
        // Incrementing counter and removing the row
        ++ans.at(i, 1u);
        x0.erase(x0.begin() + j);
        
      } else j++;
    }
  }
  
  return ans;
}

// [[Rcpp::export(rng=false)]]
IntegerVector fast_table_using_labels(
    const IntegerVector & x,
    const IntegerVector & ids
) {
  
  std::vector< int > x0(x.size());
  for (unsigned int i = 0u; i < x0.size(); ++i)
    x0[i] = x.at(i);
  IntegerVector ans(ids.size());
  
  unsigned int N = ids.size(), i, j;
  
  for (i = 0u; i < N; i++) {
    
    // Filling the first row of the output
    ans.at(i) = 0;
    
    // Looping through xsize
    j = 0u;
    while (j < x0.size()) {
      
      // If id of i is in x0, then add it!
      if (x0[j] == ids.at(i)) {
        
        // Incrementing counter and removing the row
        // that we just counted
        ++ans.at(i);
        x0.erase(x0.begin() + j);
        
      } else ++j;
    }
  }
  
  return ans;
}

typedef std::vector< std::vector<int> > stdintvec;

// [[Rcpp::export(name = ".list_offspring", rng=false)]]
ListOf<IntegerVector> list_offspring(IntegerMatrix E, int n) {
  stdintvec ans(n);
  
  for (int i = 0; i < E.nrow(); i++)
    ans.at(E.at(i, 0) - 1).push_back(E.at(i, 1));
  
  List O(n);
  for (int i = 0; i < n; i++)
    O.at(i) = Rcpp::wrap(ans.at(i));
  
  return wrap(O);
}


