#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export(name = "approx_geodesic.", rng=false)]]
arma::imat approx_geodesic(
    const arma::umat & edges,
    unsigned int nsteps = 5e3,
    bool undirected = true,
    bool warn = false
) {
  
  // Size of the tree
  int n = 1 + (int) arma::max(arma::max(edges)) ;
  
  // It is undirected
  arma::umat edges2 = edges;
  
  if (undirected) 
    edges2.insert_rows(
      edges.n_rows,
      arma::join_rows(edges.col(1u), edges.col(0u))
    );
  
  arma::colvec values(edges2.n_rows, arma::fill::ones);
  
  // Filling the matrix
  arma::imat G(n, n, arma::fill::zeros);
  
  // Filling the matrix
  for (unsigned int i = 0; i < edges2.n_rows; ++i)
    G.at(edges2.at(i, 0u), edges2.at(i, 1u)) = 1;
  
  arma::imat ans(n,n);
  ans.fill(-1);
  
  // Going through the steps
  arma::imat pG = G;
  nsteps = ((unsigned int) n) > nsteps ? nsteps : (int) n;
  for (int iter = 0; iter < (int) nsteps; ++iter) {
    
    // Computing nsteps
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        if (i != j && pG.at(i, j) != 0 && (ans.at(i, j) == -1)) 
          ans.at(i, j) = iter + 1;
      }
      
      
    // Graph power
    pG *= G;
  }
  
  // Filling diagonal with zeros
  ans.diag().zeros();
  
  return ans;
}

//' Matrix of states
//' 
//' @param P Integer scalar. Number of functions.
//' @return A matrix of size 2^P by P with all the possible
//' (0,1) combinations of functions.
//' @examples
//' states(3)
//' @export
// [[Rcpp::export(rng=false)]]
arma::imat states(
    int P
) {
  
  // Creating output matrix
  int nstates = pow(2, P);
  arma::imat ans(nstates, P);
  
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
arma::mat prob_mat(
    const arma::vec & pr
) {
  arma::mat ans(2,2);
  
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++) 
      ans.at(i,j) = 
        !i?
        (  j? pr.at(0) : (1-pr.at(0)) ):
      ( !j? pr.at(1) : (1-pr.at(1)) );
  
  return ans;
}

// [[Rcpp::export(rng = false)]]
arma::vec root_node_prob(
    double Pi,
    const arma::imat & S
) {
  // Obtaining relevant constants
  int P       = S.n_cols;
  int nstates = S.n_rows;
  
  arma::vec ans(nstates);
  ans.ones();
  
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
    const arma::ivec & pseq,
    const arma::mat  & A,
    const List & offspring
) {
  

  int P = A.n_cols;
  std::vector< int > newpseq;
  std::vector< bool > included(A.n_rows);
  

  typedef arma::ivec::const_iterator iviter;
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
