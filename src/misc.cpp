#include <Rcpp.h>
using namespace Rcpp;

//' Matrix of states
//' 
//' @param P Integer scalar. Number of functions.
//' @return A matrix of size 2^P by P with all the possible
//' (0,1) combinations of functions.
//' @examples
//' states(3)
//' @export
// [[Rcpp::export(rng=false)]]
IntegerMatrix states(
    int P
) {
  
  // Creating output matrix
  int nstates = pow(2, P);
  IntegerMatrix ans(nstates, P);
  
  // Go through states
  for (int i=0; i<nstates; i++) {
    int x=i;
    
    // Go through functions
    for (int p=0; p<P; p++) {
      ans.at(i, p) = x%2;
      x /= 2;
    }
  }
  
  return ans;
}

// Probability Matrix
// 
// Generates a 2x2 matrix with probabilities for states 0/1, with rows and columns
// corresponding to states (0,1) of parent and offspring respectively.
// 
// @param pr A numeric vector of length 2.
// 
// @return a 2x2 matrix with 0/1 probabilities, with rows and columns
// corresponding to states (0,1) of parent and offspring respectively.
// @export
// [[Rcpp::export(rng=false)]]
NumericMatrix prob_mat(
    const NumericVector & pr
) {
  NumericMatrix ans(2,2);
  
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++) 
      ans.at(i,j) = 
        !i?
        (  j? pr.at(0) : (1-pr.at(0)) ):
      ( !j? pr.at(1) : (1-pr.at(1)) );
  
  return ans;
}

// [[Rcpp::export(rng = false)]]
NumericVector root_node_prob(
    double Pi,
    const IntegerMatrix & S
) {
  // Obtaining relevant constants
  int P       = S.ncol();
  int nstates = S.nrow();
  
  NumericVector ans(nstates, 1.0);

  for (int s=0; s<nstates; s++)
    for (int p=0; p<P; p++)
      ans.at(s) *= (S.at(s,p) == 0)? (1.0 - Pi) : Pi;
  
  return ans;
}

//' Reduces the peeling sequence so that only nodes that have something to contribute
//' are included in the sequence.
//' @noRd
// [[Rcpp::export(rng=false)]]
IntegerVector reduce_pseq(
    const IntegerVector & pseq,
    const NumericMatrix  & A,
    const List & offspring
) {
  

  int P = A.ncol();
  std::vector< int > newpseq;
  std::vector< bool > included(A.nrow());
  

  typedef IntegerVector::const_iterator iviter;
  typedef IntegerVector::const_iterator Riviter;
  IntegerVector O;
  for (iviter i = pseq.begin(); i != pseq.end(); ++i) {
    
    // Default to false
    included.at(*i - 1u) = false;
    
    // First check its offspring
    if (Rf_length(offspring.at(*i - 1u))) {
      O = offspring.at(*i - 1u);
      for (Riviter o = O.begin(); o != O.end(); ++o) {
        
        // Was any offspring included?
        if (included.at(*o - 1u)) {
          included.at(*i - 1u)  = true;
          newpseq.push_back(*i);
          break;
        }
        
      }
      
      // Did we included it? If yes, then go to the next dude
      if (included.at(*i - 1u))
        continue;
        
    }
    
    // We'll try to include him looking at the annotations...
    for (int p = 0; p < P; p++) 
      if (A.at(*i - 1u, p) != 9) {
        included.at(*i - 1u) = true;
        newpseq.push_back(*i);
        break;
      }
      
  }
    
  return Rcpp::wrap(newpseq);
}
