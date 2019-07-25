#include <Rcpp.h>
#include "pruner.hpp"
using namespace Rcpp;

// This function creates pre-filled arrays
template <class T>
inline std::vector< std::vector< T > > new_vector_array(
    uint n, uint m, const T val
) {
  
  typedef std::vector< std::vector< T > > vvt;
  std::vector< T > base(m, val);
  vvt ans(n, base);
  
  return ans;
}

// Generates matrix of possible sets f states
inline std::vector< std::vector< unsigned int > > states_mat(
    uint P
  ) {
  
  // Creating output matrix
  uint nstates = (uint) pow(2, P);
  std::vector< std::vector< unsigned int > > ans = new_vector_array(nstates, P, 0u);
  
  // Go through states
  for (uint i = 0u; i < nstates; ++i) {
    uint x = i;
    
    // Go through functions
    for (uint p = 0u; p < P; ++p) {
      ans[i][p] = x%2; 
      x /= 2;
    } 
  }
  
  return ans;
}

// Replaces values of a transition matrix
inline void transition_mat(
    const std::vector< double > & pr,
    std::vector< std::vector< double > > & ans
  ) {
  
  for (uint i = 0u; i < 2u; ++i)
    for (uint j = 0u; j < 2u; ++j)
      ans[i][j] =
        !i?
        (  j? pr[0] : (1.0-pr[0]) ):
      ( !j? pr[1] : (1.0-pr[1]) );
  
  return;
  
}

// Initializes a transition matrix (2 x 2)
inline std::vector< std::vector< double > > transition_mat(
    const std::vector< double > & pr
) {
  
  std::vector< std::vector< double > > ans = new_vector_array(2u, 2u, 0.0);
  
  transition_mat(pr, ans);
  
  return ans;
}

// Computes the vector of rootnode probs depending on the set of possible
// states
inline void root_node_pr(
    std::vector< double > & Pr_root,
    double pi,
    const std::vector< std::vector< unsigned int > > & S
) {
  
  for (uint s = 0u; s < S.size(); ++s) {
    Pr_root[s] = 1.0;
    for (uint p = 0u; p < S[0].size(); ++p)
      Pr_root[s] *= (S[s][p] == 0u)? (1.0 - pi) : pi;
  }
  
  return;
  
}

/*******************************************************************************
Definition of tree data 
*******************************************************************************/

class pruner::TreeData {
  
public:
  
  uint nstates;
  uint n;
  uint nfuns;
  uint nannotated;
  double prop_type_d;
  
  // Annotations
  const pruner::vv_uint A;
  const pruner::v_uint types;
  
  // Temporal storage ----------------------------------------------------------
  vv_uint states;
  vv_dbl Pr;
  double ll;
  
  // Model parameters
  vv_dbl PSI,
    // Duplication and Speciation mu
    MU_d, MU_s;
  std::vector< vv_dbl* > MU;
  v_dbl eta, Pi;  
  
  void set_mu_d(const v_dbl & mu_d_) {transition_mat(mu_d_, this->MU_d);return;}
  void set_mu_s(const v_dbl & mu_s_) {transition_mat(mu_s_, this->MU_s);return;}
  void set_psi(const v_dbl & psi_) {transition_mat(psi_, this->PSI);return;}
  void set_eta(const v_dbl & eta_) {this->eta = eta_;return;}
  void  set_pi(double pi_) {root_node_pr(this->Pi, pi_, states);return;}
  
  // Destructor and constructor ------------------------------------------------
  ~TreeData() {};
  TreeData(const vv_uint A_, const v_uint Ntype_, uint nannotated) : A(A_), types(Ntype_) {
    
    // Initializing data
    // this->A       = A;
    // this->types   = types;
    
    // Getting meta info, and initializing containers
    this->nfuns      = A[0].size();
    this->n          = A.size();
    this->nannotated = nannotated;
    this->states     = states_mat(this->nfuns);
    this->nstates    = this->states.size();
    this->Pr         = new_vector_array(this->n, this->nstates, 1.0);
    
    // Initializing parameter containers
    eta.resize(2u, 0.0);
    Pi.resize(nstates, 0.0);
    
    MU_d.resize(2u);
    MU_d[0].resize(2u);
    MU_d[1].resize(2u);
    
    MU_s.resize(2u);
    MU_s[0].resize(2u);
    MU_s[1].resize(2u);
    
    MU.resize(2u);
    
    PSI.resize(2u);
    PSI[0].resize(2u);
    PSI[1].resize(2u);
    
    // Counting the proportion of type 0
    double increments = 1.0/this->n;
    this->prop_type_d = 0.0;
    for (auto iter = types.begin(); iter != types.end(); ++iter) {
      
      if (*iter == 0u) {
        this->prop_type_d =+ increments;
      } else if (*iter != 1u)
        stop("Values in the type of node should be either 0 or 1.");
      
    }
    
    ll = 0.0;
    
  };
};
