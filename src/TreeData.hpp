#include <Rcpp.h>
#include "pruner.hpp"
using namespace Rcpp;

class pruner::TreeData {
  
public:
  
  uint nstates;
  uint n;
  uint nfuns;
  
  // Annotations
  pruner::vv_uint A;
  pruner::v_uint Ntype;
  
  // Temporal storage ----------------------------------------------------------
  vv_uint states;
  vv_dbl PSI, MU;
  vv_dbl Pr;
  v_dbl Pi;
  
  // Model parameters
  v_dbl mu, psi, eta;
  double pi;
  double ll;
  
  void  set_mu(const v_dbl & mu) {prob_mat(mu, this->MU);return;}
  void set_psi(const v_dbl & psi) {prob_mat(psi, this->PSI);return;}
  void set_eta(const v_dbl & eta) {this->eta = eta;return;}
  void  set_pi(double pi) {root_node_pr(this->Pi, pi, states);return;}
  
  // Destructor and constructor ------------------------------------------------
  ~TreeData() {};
  TreeData(vv_uint A, v_uint Ntype) {
    
    // Initializing data
    this->A       = A;
    this->Ntype   = Ntype;
    
    // Getting meta info, and initializing containers
    this->nfuns   = A[0].size();
    this->n       = A.size();
    this->states  = states_mat(this->nfuns);
    this->nstates = this->states.size();
    this->Pr      = new_vector_array(this->n, this->nstates, 1.0);
    
    // Initializing parameter containers
    mu.resize(2u, 0.0);
    psi.resize(2u, 0.0);
    eta.resize(2u, 0.0);
    Pi.resize(nstates, 0.0);
    pi = 0.0;
    
    MU  = prob_mat(mu);
    PSI = prob_mat(psi);
    
    ll = 0.0;
    
  };
};
