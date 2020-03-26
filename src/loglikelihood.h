#include "pruner.hpp"
#include "TreeData.hpp" // TreeData definition

#ifndef APHYLO_LOGLIKELIHOOD_H
#define APHYLO_LOGLIKELIHOOD_H 1

void likelihood(
    pruner::sptr_treedata D,
    pruner::TreeIterator & n
) {
  
#ifdef DEBUG_LIKELIHOOD
  printf("Entering likelihood at node %d with pruneseq:\n", *n);
  pruner::v_uint xx = n.tree->get_postorder();
  for (auto iter = xx.begin(); iter != xx.end(); ++iter)
    printf("%i, ", *iter);
  printf("\n");
#endif
  
  if (n.is_tip()) {
    
    // Iterating through the states
    pruner::uint s, p;
    for (s = 0u; s < D->states.size(); ++s) {
      
      // Throught the functions
      D->Pr[*n][s] = 1.0; // Initializing
      for (p = 0u; p < D->nfuns; ++p) {
        
        // ETA PARAMETER
        if (D->A[*n][p] == 9u && (D->eta[0u] >= 0.0)) {
          
          D->Pr[*n][s] *=
            (1.0 - D->eta[0u]) * D->PSI[D->states[s][p]][0u] +
            (1.0 - D->eta[1u]) * D->PSI[D->states[s][p]][1u]
          ;
          
        } else {
          
          // Unnanotated leafs should be skipped in this situation. This mostly
          // happens if we are not using the eta parameter and during cv. At
          // that time some annotations are dropped without modifying the 
          // pruning sequence.
          if (D->A[*n][p] == 9u)
            continue;
          
          if (D->eta[0u] >= 0.0) {
            
            D->Pr[*n][s] *= D->PSI[D->states[s][p]][D->A[*n][p]]*
              D->eta[D->A[*n][p]];
            
          } else {
            
            D->Pr[*n][s] *= D->PSI[D->states[s][p]][D->A[*n][p]];
            
          }
          
        }
        
      }
      
    }
    
  } else {
    
    D->MU[0] = &(D->MU_d);
    D->MU[1] = &(D->MU_s);
    
    std::vector< unsigned int >::const_iterator o_n;
    pruner::uint s_n, p_n, s;
    double offspring_ll, s_n_sum;
    
    // Looping through states
    for (s = 0u; s < D->nstates; ++s) {
      
      // Now through offspring
      D->Pr[*n][s] = 1.0;
      for (o_n = n.begin_off(); o_n != n.end_off(); ++o_n) {
        
        // Offspring state integration
        offspring_ll = 0.0;
        for (s_n = 0u; s_n < D->nstates; ++s_n) {
          
          s_n_sum = 1.0;
          for (p_n = 0u; p_n < D->nfuns; ++p_n)
            // s_n_sum *= (D->MU[D->types[*n]]).at(D->states[s][p_n]).at(D->states[s_n][p_n]);
            s_n_sum *= (D->types[*n] == 0u)?
          D->MU_d[D->states[s][p_n]][D->states[s_n][p_n]] :
            D->MU_s[D->states[s][p_n]][D->states[s_n][p_n]];
          
          // Multiplying by off's probability
          offspring_ll += (s_n_sum) * D->Pr[*o_n][s_n];
          
        }
        
        // Getting the joint conditional.
        D->Pr[*n][s] *= offspring_ll;
        
      }
      
    }
    
    // Computing the joint likelihood
    if (*n == n.back()) {
      D->ll = 0.0;
      for (s = 0; s < D->nstates; ++s) 
        D->ll += D->Pi[s] * D->Pr[*n][s];
      D->ll = log(D->ll);
    }
    
  }
  
  
  return;
  
}

#endif

