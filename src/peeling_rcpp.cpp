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

//' Leaf probabilities
//' 
//' @param Z A matrix of size \eqn{N\times P}{N*P} with values \code{(0,1,9)}.
//' @param S A matrix of size \eqn{2^P*P} as returned by \code{\link{states}}.
//' @param psi A numeric vector of length 2.
//' @param noffsprings A numeric vector of length \eqn{N}. Number of offsprings
//' per node.
//' 
//' @export
// [[Rcpp::export]]
arma::mat leaf_prob(
    const arma::imat & Z,
    const arma::imat & S,
    const arma::vec  & psi,
    const arma::ivec & noffsprings
  ) {
  
  // Obtaining relevant constants
  int P       = S.n_cols;
  int N       = Z.n_rows;
  int nstates = S.n_rows;
  
  // Checking dimmensions
  if (P != (int) Z.n_cols)
    stop("Z must have the same number of columns than S.");
  
  if (2u != psi.n_elem)
    stop("psi must be of length 2.");
  
  // Creating the output matrix (probabilities)
  arma::mat ans(N,nstates);
  ans.ones();

  for (int i=0; i<N; i++)
    // Only compute for offsprings (and nodes with no NA)
    if (!noffsprings.at(i)) 
      for (int s=0; s<nstates; s++)
        for (int p=0; p<P; p++) {
          
          // If missing (no experimental data)
          if (Z.at(i,p) == 9)
            continue;
          
          ans.at(i, s) *= (Z.at(i, p) == S.at(s, p))? (1 - psi.at(0)): psi.at(1);
          
        }
      
  return ans;
}

//' Gain/Loss probabilities
//' 
//' Generates a 2x2 matrix with gain loss probabilities, with rows and columns
//' corresponding to states (0,1) of parent and offspring respectively.
//' 
//' @param mu A numeric vector of length 2. Gain and Loss probailities.
//' 
//' @return a 2x2 matrix with gain loss probabilities, with rows and columns
//' corresponding to states (0,1) of parent and offspring respectively.
//' @export
// [[Rcpp::export]]
arma::mat gain_loss_prob(
    const arma::vec & mu
  ) {
  arma::mat ans(2,2);
  
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++) 
      ans.at(i,j) = !i ? ( j? mu.at(0) : (1-mu.at(0)) ) :
      ( !j? mu.at(1) : (1-mu.at(1)) );
    
  return ans;
}

//' Root node probabilities
//' 
//' Generates a vector of length 2^P with the root node probabilities
//' per state
//' @export
// [[Rcpp::export]]
arma::vec root_node_prob(
    const arma::vec  & pi,
    const arma::imat & S
) {
  // Obtaining relevant constants
  int P       = S.n_cols;
  int nstates = S.n_rows;
  
  arma::vec ans(nstates);
  ans.ones();
  
  for (int s=0; s<nstates; s++)
    for (int p=0; p<P; p++)
      ans.at(s) *= pi.at(S.at(s,p));
  
  return ans;
}

//' Internal node probabilities
//' 
//' @param Pr Probabilities (already with leaf probs).
//' @param M Gain/Loss probabilities (see equation 4 of math.pdf)
//' @param S States
//' @param pi Root node state probabilities (as returned by root_node_prob)
//' @param noffsprings Number of offsprings
//' @param offsprings List of offsprings
//' @export
// [[Rcpp::export]]
arma::mat internal_prob(
  arma::mat          Pr,
  const arma::mat  & M,
  const arma::imat & S,
  const arma::ivec & noffsprings,
  const List       & offsprings
) {
  
  // Obtaining relevant constants
  int P       = S.n_cols;
  int N       = Pr.n_rows;
  int nstates = S.n_rows;
  
  for (int n=(N-1); n>=0; n--) {
    
    // Only for root nodes
    if (!noffsprings.at(n))
      continue;
    
    // Parent node states integration
    for (int s=0; s<nstates; s++) {
    
      // Obtaining list of offsprings
      IntegerVector O(offsprings.at(n));
      
      // Loop through offsprings
      double offsprings_joint_likelihood = 1.0;
      for (int o_n=0; o_n<noffsprings.at(n) ; o_n++) {

        // Offspring states integration
        double offspring_likelihood = 0.0;
        for (int s_n=0; s_n<nstates; s_n++) {
          double s_n_sum = 1.0;
          
          // Loop through functions
          for (int p=0; p<P; p++) 
            s_n_sum *= M.at(s, s_n);
          
          // Multiplying by prob of offspring
          offspring_likelihood += (s_n_sum *Pr.at(O.at(o_n), s_n));
          
        }
        
        // Multiplying with other offsprings
        offsprings_joint_likelihood *= offspring_likelihood;
        
      }
      
      // Adding state probability
      Pr.at(n, s) = offsprings_joint_likelihood;
      
    }
  }
  
  return Pr;
}

//' Log likelihood
//' 
//' @export
// [[Rcpp::export]]
List LogLike(
    const arma::imat & Z,
    const List       & offsprings,
    const arma::ivec & noffsprings,
    const arma::vec  & psi,
    const arma::vec  & mu,
    const arma::vec  & Pi
) {

  // Obtaining States, PSI, Gain/Loss probs, and root node probs
  arma::imat S  = states(Z.n_cols);
  int nstates   = (int) S.n_cols;
  arma::mat PSI = leaf_prob(Z, S, psi, noffsprings);
  arma::mat M   = gain_loss_prob(mu);
  arma::vec PiP = root_node_prob(Pi, S);
  
  // Computing likelihood
  arma::mat Pr  = internal_prob(PSI, M, S, noffsprings, offsprings);
  
  // We only use the root node
  double ll = 0.0;
  for (int s = 0; s<nstates; s++)
    ll += log(PiP.at(s)*Pr.at(0, s));
  
  // return ll;
  return List::create(
    _["S"]   = S,
    _["PI"]  = PiP,
    _["PSI"] = PSI,
    _["Pr"]   = Pr,
    _["ll"]  = ll
  );
  
}

