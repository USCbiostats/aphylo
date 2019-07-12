#include "typedefs.hpp"

// This function creates pre-filled arrays
template <class T>
inline std::vector< std::vector< T > > new_vector_array(uint n, uint m, const T val) {
  
  typedef std::vector< std::vector< T > > vvt;
  std::vector< T > base(m, val);
  vvt ans(n, base);
  
  return ans;
}

// Takes a vec<vec> and resets it to a default value
template <class T>
inline void reset_values(std::vector< std::vector< T > > & vec, const T val) {
  
  for (auto i = vec.begin(); i != vec.end(); ++i)
    std::fill(i->begin(), i->end(), val);
  
  return;
  
}

inline std::vector< std::vector< unsigned int > > states_mat(uint P) {
  
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

inline void prob_mat(
    const std::vector< double > & pr,
    std::vector< std::vector< double > > & ans) {
  
  for (uint i = 0u; i < 2u; ++i)
    for (uint j = 0u; j < 2u; ++j)
      ans[i][j] =
        !i?
        (  j? pr[0] : (1.0-pr[0]) ):
      ( !j? pr[1] : (1.0-pr[1]) );
  
  return;
  
}

inline std::vector< std::vector< double > > prob_mat(const std::vector< double > & pr) {
  
  std::vector< std::vector< double > > ans = new_vector_array(2u, 2u, 0.0);
  
  prob_mat(pr, ans);
  
  return ans;
}

inline std::vector< std::vector< double > > probabilities(uint N, uint nstates) {
  
  return new_vector_array(N, nstates, 1.0);
  
}

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
