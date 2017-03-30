// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
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
      ans.at(i,j) = !i ? ( j? pr.at(0) : (1-pr.at(0)) ) :
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
//' @templateVar Z 1
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
    const arma::imat & Z,
    const arma::vec  & mu,
    const arma::vec  & psi,
    const arma::imat & S,
    const arma::ivec & noffspring,
    const List       & offspring
) {
  
  // Obtaining relevant constants
  int P         = S.n_cols;
  int N         = Z.n_rows;
  int nstates   = S.n_rows;
  arma::mat M   = prob_mat(mu);
  arma::mat PSI = prob_mat(psi);
  
  arma::mat Pr(N, nstates);
  Pr.ones();
  
  for (int n=(N-1); n>=0; n--) {
    
    // Only for internal nodes
    if (!noffspring.at(n)) {
      for (int s=0; s<nstates; s++)
        for (int p=0; p<P; p++) {
          
          // If missing (no experimental data)
          if (Z.at(n,p) == 9)
            continue;
          
          Pr.at(n, s) *= PSI.at(S.at(s, p), Z.at(n, p));
          
        }
        continue;
    }
      
    
    // Parent node states integration
    for (int s=0; s<nstates; s++) {
      
      // Obtaining list of offspring <- this can be improved (speed)
      // can create an std vector of size n
      IntegerVector O(offspring.at(n));
      
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
//' a particular dataset. The arguments \code{Z}, \code{offspring}, and
//' \code{noffspring} should be as those returned by \code{\link{get_offspring}}.
//' For complete Maximum Likelihood Estimation see \code{\link{mle}}.
//' 
//' @template parameters
//' @templateVar Z 1
//' @templateVar offspring 1
//' @templateVar noffspring 1
//' @templateVar psi 1
//' @templateVar mu 1
//' @templateVar Pi 1
//' @param verb_ans Logical scalar. When \code{FALSE} (default) the function
//' returns a list with a single scalar (the log-likelihood).
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
    const arma::imat & Z,
    const List       & offspring,
    const arma::ivec & noffspring,
    const arma::vec  & psi,
    const arma::vec  & mu,
    const arma::vec  & Pi,
    bool verb_ans = false
) {

  // Obtaining States, PSI, Gain/Loss probs, and root node probs
  arma::imat S  = states(Z.n_cols);
  int nstates   = (int) S.n_rows;

  // Computing likelihood
  arma::mat Pr  = probabilities(Z, mu, psi, S, noffspring, offspring);

  // We only use the root node
  double ll = 0.0;

  arma::vec PiP = root_node_prob(Pi, S);  
  for (int s = 0; s<nstates; s++)
    ll += log(PiP.at(s)*Pr.at(0, s));
  
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

