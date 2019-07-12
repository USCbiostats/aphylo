#include <Rcpp.h>
#include "pruner.hpp"
#include "likelihood.hpp"
#include "TreeData.hpp" // TreeData definition
using namespace Rcpp;

void likelihood(
    pruner::sptr_treedata D,
    pruner::TreeIterator & n
) {
  
  if (n.is_leaf()) {
    
    // Iterating through the states
    for (uint s = 0u; s < D->states.size(); ++s) {
      
      // Throught the functions
      D->Pr[*n][s] = 1.0; // Initializing
      for (uint p = 0u; p < D->nfuns; ++p) {
        
        // ETA PARAMETER
        if (D->A[*n][p] == 9u) {
          
          D->Pr[*n][s] *=
            (1.0 - D->eta[0u]) * D->PSI[D->states[s][p]][0u] +
            (1.0 - D->eta[1u]) * D->PSI[D->states[s][p]][1u]
            ;
          
        } else {
          D->Pr[*n][s] *= D->PSI[D->states[s][p]][D->A[*n][p]]*
            D->eta[D->A[*n][p]];
        }
          
      }
      
    }
    
  } else {
    
    std::vector< unsigned int >::const_iterator o_n;
    uint p, s_n, p_n;
    double offspring_ll, s_n_sum;
    
    // Looping through states
    for (uint s = 0u; s < D->nstates; ++s) {
      
      // Now through offspring
      D->Pr[*n][s] = 1.0;
      for (o_n = n.begin_off(); o_n != n.end_off(); ++o_n) {
        
        // Offspring state integration
        offspring_ll = 0.0;
        for (s_n = 0u; s_n < D->nstates; ++s_n) {
          
          s_n_sum = 1.0;
          for (p_n = 0u; p_n < D->nfuns; ++p_n)
            s_n_sum *= D->MU[D->states[s][p]][D->states[s_n][p]];
          
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
      for (uint s = 0; s < D->nstates; ++s) 
        D->ll += D->Pi[s] * D->Pr[*n][s];
      D->ll = log(D->ll);
    }
    
  }
  
  
  return;
  
}

// Tree constructor ------------------------------------------------------------

// [[Rcpp::export(rng = false)]]
SEXP Tree_new(
    const std::vector< std::vector< unsigned int > > & edgelist,
    const std::vector< std::vector< unsigned int > > & A,
    const std::vector< unsigned int >  & Ntype
) {
  
  // Initializing the tree
  uint res;
  Rcpp::XPtr< pruner::Tree > xptr(new pruner::Tree(edgelist[0], edgelist[1], res), true);
  
  xptr->args = std::make_shared< pruner::TreeData >(A, Ntype);
  xptr->fun  = likelihood;
  
  return wrap(xptr);
}

// Methods ---------------------------------------------------------------------

// [[Rcpp::export(name = ".LogLike2", rng = false)]]
List LogLike2(
    SEXP tree_ptr,
    const std::vector< double > & mu,
    const std::vector< double > & psi,
    const std::vector< double > & eta,
    const double & pi,
    bool verb = true
) {
  Rcpp::XPtr< pruner::Tree > p(tree_ptr);
  
  // Setting the parameters
  p->args->set_mu(mu);
  p->args->set_psi(psi);
  p->args->set_pi(pi);
  p->args->set_eta(eta);
  
  p->prune_postorder();
  
  if (verb)
    return List::create(
      _["Pr"] = wrap(p->args->Pr), _["ll"] = wrap(p->args->ll)
    );
  else
    return List::create(_["ll"] = wrap(p->args->ll));
}

// [[Rcpp::export]]
std::vector< std::vector< unsigned int > > Tree_get_offspring(const SEXP & tree_ptr) {
  
  Rcpp::XPtr< pruner::Tree > p(tree_ptr);
  return p->get_offspring();
  
}

// [[Rcpp::export]]
std::vector< std::vector< unsigned int > > Tree_get_parents(const SEXP & tree_ptr) {
  
  Rcpp::XPtr< pruner::Tree > p(tree_ptr);
  return p->get_parents();
  
}


/***R
set.seed(1)
dat <- aphylo::raphylo(50)

A <- rbind(dat$tip.annotation, dat$node.annotation)[,1]

mu  <- c(.1, .05)
psi <- c(.2, .07)
eta <- c(.8, .9)
Pi  <- .5

tree_ptr <- aphylo:::new_Tree(
  edgelist = with(dat$tree, list(edge[,1] - 1, edge[,2] - 1)),
  A        = as.list(A),
  Ntype    = A
)

aphylo_ll <- aphylo::LogLike

microbenchmark::microbenchmark(
  new = aphylo:::.LogLike2(tree_ptr, mu = mu, psi = psi, eta = eta, pi = Pi),
  old = aphylo::LogLike(dat, psi = psi, mu = mu, Pi = Pi, eta = eta)$ll,
  times = 100, unit = "relative"
)


*/
