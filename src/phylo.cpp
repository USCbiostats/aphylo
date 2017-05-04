// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "misc.h"
using namespace Rcpp;

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

// Root node probabilities
// 
// @template parameters
// @templateVar Pi 1
// @templateVar S 1
// 
// Generates a vector of length 2^P with the root node probabilities
// per state
// @export
// [[Rcpp::export]]
arma::vec root_node_prob(
    const arma::vec  & Pi,
    const arma::imat & S
) {
  // Obtaining relevant constants
  int P       = S.n_cols;
  int nstates = S.n_rows;
  
  arma::vec ans(nstates);
  ans.ones();
  
  for (int s=0; s<nstates; s++)
    for (int p=0; p<P; p++)
      ans.at(s) *= Pi.at(S.at(s,p));
  
  return ans;
}

//' State probabilities
//' 
//' Compute the state probabilities for each node in the tree using the peeling
//' algorithm. This function is the horse-power of the function \code{\link{LogLike}}
//' so it is not intended to be used directly.
//' 
//' @template parameters
//' @templateVar annotations 1
//' @templateVar mu 1
//' @templateVar psi 1
//' @templateVar S 1
//' @templateVar noffspring 1
//' @templateVar offspring 1
//' 
//' @export
//' @return A numeric matrix of size \eqn{n\times 2^P}{n * 2^P} with state
//' probabilities for each node.
//' 
// [[Rcpp::export]]
arma::mat probabilities(
    const arma::imat & annotations,
    const arma::vec  & mu,
    const arma::vec  & psi,
    const arma::imat & S,
    const arma::ivec & noffspring,
    const List       & offspring
) {
  
  // Obtaining relevant constants
  int P         = S.n_cols;
  int N         = annotations.n_rows;
  int nstates   = S.n_rows;
  arma::mat M   = prob_mat(mu);
  arma::mat PSI = prob_mat(psi);
  
  arma::mat Pr(N, nstates);
  Pr.ones();
  
  for (int n=(N-1); n>=0; n--) {
    // Rprintf("Looping in n=%i\n", n);
    // Only for internal nodes
    if (!noffspring.at(n)) {
      for (int s=0; s<nstates; s++)
        for (int p=0; p<P; p++) {
          
          // If missing (no experimental data)
          if (annotations.at(n,p) == 9)
            continue;
          
          Pr.at(n, s) *= PSI.at(S.at(s, p), annotations.at(n, p));
          
        }
        continue;
    }
      
    // Obtaining list of offspring <- this can be improved (speed)
    // can create an std vector of size n
    IntegerVector O(offspring.at(n));
    
    // Parent node states integration
    for (int s=0; s<nstates; s++) {
      
      // Loop through offspring
      double offspring_joint_likelihood = 1.0;
      for (int o_n=0; o_n<noffspring.at(n) ; o_n++) {
        
        // Offspring states integration
        double offspring_likelihood = 0.0;
        for (int s_n=0; s_n<nstates; s_n++) {
          double s_n_sum = 1.0;
          
          // Loop through functions
          for (int p=0; p<P; p++) 
            s_n_sum *= M.at(S.at(s,p), S.at(s_n,p));
          
          // Multiplying by prob of offspring
          offspring_likelihood += ( s_n_sum * Pr.at(O.at(o_n), s_n) );
          
        }
        
        // Multiplying with other offspring
        offspring_joint_likelihood *= offspring_likelihood;
        
      }
      
      // Adding state probability
      Pr.at(n, s) = offspring_joint_likelihood;
      
    }
  }
  
  return Pr;
}

//' Computes Log-likelihood
//' 
//' This function computes the log-likelihood of the chosen parameters given
//' a particular dataset. The arguments \code{annotations}, \code{offspring}, and
//' \code{noffspring} should be as those returned by \code{\link{new_aphylo}}.
//' For complete Maximum Likelihood Estimation see \code{\link{mle}}.
//' 
//' @template parameters
//' @templateVar annotations 1
//' @templateVar offspring 1
//' @templateVar noffspring 1
//' @templateVar psi 1
//' @templateVar mu 1
//' @templateVar Pi 1
//' @param verb_ans Logical scalar. When \code{FALSE} (default) the function
//' returns a list with a single scalar (the log-likelihood).
//' @param check_dims Logical scalar. When \code{TRUE} (default) the function
//' checks the dimmension of the passed parameters.
//' 
//' @details
//' The parameters to estimate are described as follows:
//' \enumerate{
//' \item{\code{psi}: A vector of length 2 with \eqn{\psi_0}{psi[0]} and
//' \eqn{\psi_1}{psi[1]}, which are the misclassification probabilities fo
//' \eqn{s_p=0}{s[p]=0} and \eqn{s_p=1}{s[p]=1}
//' respectively.}
//' \item{\code{mu}: A vector of length 2 with \eqn{\mu_0}{mu[0]} and
//' \eqn{\mu_1}{mu[1]} which are the gain and loss probabilities respectively.}
//' \item{\code{Pi}: A vector of length 2 with \eqn{\pi_0}{pi[0]}} and
//' \eqn{\pi_1}{pi[1]}, which for now is specified as \eqn{(1 - \pi_1)}{(1 - pi[0]),
//' which holds the root node probabilities.}
//' }
//' @return A list of class \code{phylo_LogLik} with the following elements:
//' \item{S}{An integer matrix of size \eqn{2^p\times p}{2^p * p} as returned
//' by \code{\link{states}}.}
//' \item{Pr}{A numeric matrix of size \eqn{G\times 2^p}{G * 2^p} with node/state
//' probabilities.}
//' \item{ll}{A numeric scalar with the log-likelihood value given the chosen
//' parameters.}
//' @export
// [[Rcpp::export]]
List LogLike(
    const arma::imat & annotations,
    const List       & offspring,
    const arma::ivec & noffspring,
    const arma::vec  & psi,
    const arma::vec  & mu,
    const arma::vec  & Pi,
    bool verb_ans = false,
    bool check_dims = true
) {

  // Checking dimmensions
  if (check_dims) {
    
    bool dims_are_ok = true;
    
    // Data dims
    int n_annotations = annotations.n_rows;
    int n_offspring   = offspring.size();
    int n_noffspring  = noffspring.size();
    
    if (n_annotations != n_offspring) {
      warning("-annotations- and -offspring- have different lengths.");
      dims_are_ok = false;
    }
      
      
    if (n_annotations != n_noffspring) {
      warning("-annotations- and -noffspring- have different lengths.");
      dims_are_ok = false;
    }
      
    
    if (n_noffspring != n_noffspring) {
      warning("-offspring- and -noffspring- have different lengths.");
      dims_are_ok = false;
    }
      
    // Parameters dims
    int n_psi = psi.size();
    int n_mu  = mu.size();
    int n_Pi  = Pi.size();
    
    if (n_psi != 2) {
      warning("-psi- must be a vector of size 2.");
      dims_are_ok = false;
    }
    
    if (n_mu != 2) {
      warning("-mu- must be a vector of size 2.");
      dims_are_ok = false;
    }
    
    if (n_Pi != 2) {
      warning("-Pi- must be a vector of size 2.");
      dims_are_ok = false;
    }
    
    // Return with error
    if (!dims_are_ok)
      stop("Check the size of the inputs.");
    
  }
  
  // Obtaining States, PSI, Gain/Loss probs, and root node probs
  arma::imat S  = states(annotations.n_cols);
  int nstates   = (int) S.n_rows;

  // Computing likelihood
  arma::mat Pr  = probabilities(annotations, mu, psi, S, noffspring, offspring);

  // We only use the root node
  double ll = 0.0;

  arma::vec PiP = root_node_prob(Pi, S);  
  for (int s = 0; s<nstates; s++)
    ll += PiP.at(s)*Pr.at(0, s);
  
  ll = std::log(ll);
  
  // return ll;
  if (verb_ans) {
    
    List ans = List::create(
      _["S"]   = S,
      _["Pr"]   = Pr,
      _["ll"]  = ll
    );
    ans.attr("class") = "phylo_LogLik";
    return(ans);
    
  } else {
    
    List ans = List::create(_["ll"] = ll);
    ans.attr("class") = "phylo_LogLik";
    return(ans);
    
  }
  
}

// [[Rcpp::export]]
double predict_fun(
  unsigned int i,
  unsigned int p,
  unsigned int di0,
  const arma::imat & annotations,
  const List       & offspring,
  const arma::ivec & noffspring,
  const arma::vec  & psi,
  const arma::vec  & mu,
  const arma::vec  & Pi
) {
  
  arma::imat annotations_filled(annotations);
  
  //----------------------------------------------------------------------------
  // Compute likelihood of a_i = 0
  annotations_filled.at(i, p) = 0;
  double likelihood_given_ai_0 = exp(
    as< double >(LogLike(annotations_filled, offspring, noffspring, psi, mu, Pi, false)["ll"])
    );
  
  // Compute likelihood of a_i = 1
  annotations_filled.at(i, p) = 1;
  double likelihood_given_ai_1 = exp(
    as< double >(LogLike(annotations_filled, offspring, noffspring, psi, mu, Pi, false)["ll"])
    );
  
  // Pr(a_i = 1 | Tree Structure only) -----------------------------------------
  arma::mat MU = prob_mat(mu);
  
  // Rasing it to the power of di0
  for (int i = 1; i < di0; i++)
    MU = MU * MU;
  
  double Pr_ai_1 = Pi.at(1) * MU.at(1, 1) + Pi.at(0) * MU.at(0, 1);
  
  // Returning:
  return likelihood_given_ai_1 / (
      likelihood_given_ai_1 + likelihood_given_ai_0*(1 - Pr_ai_1)/Pr_ai_1
  );
}

// [[Rcpp::export]]
arma::mat predict_funs(
  const arma::uvec & ids,
  const arma::umat & edges,
  const arma::imat & annotations,
  const List       & offspring,
  const arma::ivec & noffspring,
  const arma::vec  & psi,
  const arma::vec  & mu,
  const arma::vec  & Pi
) {
  
  unsigned int n = ids.size(), i;
  unsigned int P = annotations.n_cols, p;
  arma::mat ans(n, P);
  
  // Computing geodesic
  arma::umat G = approx_geodesic(edges, 1e3, true, false);
  
  for (i = 0u; i < n; i++)
    for (p = 0u; p < P; p++) {
      ans.at(i, p) = predict_fun(
        ids.at(i), p, G.at(ids.at(i), 0), annotations, offspring, noffspring, psi, mu, Pi
      );
    }
      
  
  return ans;
}



