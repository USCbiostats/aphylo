// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// This function generates a matrix of size 2^P by P with all the possible
// (0,1) combinations of functions.
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

// [[Rcpp::export]]
arma::mat leaf_prob(
    const arma::imat & Z,
    const arma::imat & S,
    const arma::vec & psi,
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

// [[Rcpp::export]]
List list_offsprings(
    const IntegerVector & id,
    const IntegerVector & parent
  ) {
  
  int N = id.length();
  std::vector< IntegerVector > ans_tmp(N);
  
  // Listing offsprings by parent
  for (int i=0; i<N; i++) {
    if (parent.at(i) == NA_INTEGER) continue;
    ans_tmp.at(parent.at(i)).push_back(i);
  }
  
  // Coercing into a List
  List ans(N);
  for (int i=0; i<N; i++)
    ans.at(i) = ans_tmp.at(i);
  
  return ans;
}

//' @param Pr Probabilities (already with leaf probs).
//' @param S States
//' @param mu Gain/Loss probabilities
//' @param noffsprings Number of offsprings
//' @param offsprings List of offsprings
// [[Rcpp::export]]
arma::mat internal_prob(
  arma::mat & Pr,
  const arma::imat & S,
  const arma::vec & mu,
  const arma::ivec & noffsprings,
  const List & offsprings
) {
  
  // Obtaining relevant constants
  int P       = S.n_cols;
  int N       = Pr.n_rows;
  int nstates = S.n_rows;
  
  for (int i=0; i<N; i++) {
    if (!noffsprings.at(i))
      continue;
    
    // Loop through the columns of Pr
    for (int s=0; s<nstates; s++) {
    
      // Obtaining list of offsprings
      Rprintf("In\n");
      Rprintf("i:%d\n", i);
      IntegerVector O(offsprings.at(i));
      Rprintf("Out\n");
      
      // Loop through offsprings
      for (int o=0; o<noffsprings.at(i) ; o++) {
        double mutatepr = 0;

        // Integrating
        for (int subs=0; subs<nstates; subs++) {
          
          // Loop through functions
          double pi = 1;
          for (int p=0; p<P; p++) {
            // The parent doesn't have the function
            if (!S.at(s, p)) {
              pi *=  S.at(subs, p) ? mu.at(0) : (1-mu.at(0));
            } else {
              pi *= !S.at(subs, p) ? mu.at(1) : (1-mu.at(1));
            }
          }
          
          // Multiplying by misspecification prob of offspring
          mutatepr += Pr.at(O.at(o), subs)*pi;
        }
        
        // Adding it up
        Pr.at(i, s) *= mutatepr; 
      }
    }
  }
  return Pr;

}

/***R


# Reading data
dag <- read.table("../data/pthr11848.dag.txt", sep="\t", header = FALSE,
                  col.names = c("NodeId", "TypeId", "ParentId"))

dat <- read.table("../data/pthr11848_sorted.txt", sep = "\t", header = FALSE)
colnames(dat) <- c(sprintf("f%02d",1:(ncol(dat)-1)), "LeafId")

# Listing offsprings
offsprings <- lapply(dat$LeafId, function(x) {
  y <- which(dag$ParentId == x)
})

# Substracting one so we canuse it in C++
offsprings <- lapply(offsprings, function(x) {
  if (length(x)) x-1
  else x
})

# Checking who is parent
noffsprings <- sapply(offsprings, length)

# Parameters
psi <- c(0.020,0.010)
mu  <- c(0.004,.001)
S   <- states(ncol(dat)-1)

# Computing leaf probabilities
ans0 <- leaf_prob(as.matrix(dat[,-4]), S, psi, noffsprings)
ans1 <- internal_prob(ans0, S, mu, noffsprings, offsprings)

*/