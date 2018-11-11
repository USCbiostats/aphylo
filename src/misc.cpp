#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//' Approximation of Geodesic distances using Matrix Powers
//' 
//' Given an adjacency matrix \eqn{A}, the geodesic can be approximated using
//' its powers, since each \eqn{(i,j)} element of \eqn{A^p} corresponds to the
//' number of \eqn{p} length steps between nodes \eqn{i} and \eqn{j}.
//' 
//' @template edges
//' @param nsteps Integer scalar. Number of maximum steps for the approximation.
//' @param warn Logical scalar. When \code{TRUE} shows a warning after no further
//' steps are needed.
//' @param undirected Logical scalar. When \code{TRUE} (default), the edgelist is treated
//' as undirected (see details).
//' @return A square matrix of size \eqn{n} with the shortest path between each
//' pair of nodes.
//' 
//' @details
//' When \code{undirected = TRUE}, the function extends \code{edges} such that
//' \code{edges = rbind(edges, edges[,2:1])}.
//' 
//' @author George G. Vega Yon
//' @references
//' This is a modified version of the function of the same name in the
//' R package \CRANpkg{netdiffuseR}.
//' @export
//[[Rcpp::export]]
arma::umat approx_geodesic(
    const arma::umat & edges,
    unsigned int nsteps = 5e3,
    bool undirected = true,
    bool warn = false
) {
  
  // Size of the tree
  int n = 1 + (int) arma::max(arma::max(edges));
  
  // It is undirected
  arma::umat edges2 = edges;
  
  if (undirected) 
    edges2.insert_rows(
      edges.n_rows,
      arma::join_rows(edges.col(1u), edges.col(0u))
    );
  
  arma::colvec values(edges2.n_rows, arma::fill::ones);
  
  // Filling the matrix
  arma::sp_mat G(edges2.t(), values, n, n, true, false);
  arma::umat ans(n,n, arma::fill::zeros);
  
  typedef arma::sp_mat::const_iterator spiter;
  
  // Going through the steps
  arma::sp_mat pG = G;
  arma::sp_mat G0 = G;
  int change_count = 0;
  nsteps++;
  for (unsigned int i=1u; i<nsteps; i++) {
    
    // Computing nsteps
    
    for (spiter it = pG.begin(); it != pG.end(); ++it)
      if (ans.at(it.row(), it.col()) == 0u)
        ans.at(it.row(), it.col()) += i,
          ++change_count;
      
      // Was there any change?
      if (!change_count) {
        if (warn)
          warning("The algorithm stopped at %i iterations.", i);
        break;
      } else change_count = 0;
      
      // Graph power
      pG *= G0;
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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

//' Reduces the peeling sequence so that only nodes that have something to contribute
//' are included in the sequence.
//' @noRd
// [[Rcpp::export]]
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