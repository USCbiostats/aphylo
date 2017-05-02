#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//' Approximation of Geodesic distances using Matrix Powers
//' 
//' Given an adjacency matrix \eqn{A}, the geodesic can be approximated using
//' its powers, since each \eqn{(i,j)} element of \eqn{A^p} corresponds to the
//' number of \eqn{p} length steps between nodes \eqn{i} and \eqn{j}.
//' 
//' @template parameters
//' @templateVar edges 1
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
    unsigned int nsteps = 1e3,
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
    
    for (spiter it = pG.begin(); it != pG.end(); it ++)
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