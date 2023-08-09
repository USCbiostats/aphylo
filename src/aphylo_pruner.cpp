#include <Rcpp.h>
#include "pruner.hpp"
#include "TreeData.hpp" // TreeData definition
#include "loglikelihood.h"
using namespace Rcpp;

// #define DEBUG_LIKELIHOOD

// Tree constructor ------------------------------------------------------------

// [[Rcpp::export(rng = false)]]
SEXP new_aphylo_pruner_cpp(
    const std::vector< std::vector< unsigned int > > & edgelist,
    const std::vector< std::vector< unsigned int > > & A,
    const std::vector< unsigned int >  & types,
    unsigned int nannotated
) {
  
  // Initializing the tree
  pruner::uint res;
  Rcpp::XPtr< AphyloPruner > xptr(
      new AphyloPruner(A, types, nannotated, edgelist[0], edgelist[1], res),
      true);
  
  if (res != 0u)
    stop(
      "An error of code %d happened while creating the pruner::Tree object.",
      res
    );
  

  xptr.attr("class") = "aphylo_pruner";
  
  return xptr;
}

// Methods ---------------------------------------------------------------------
// [[Rcpp::export]]
int sizeof_pruner(SEXP ptr) {
  Rcpp::XPtr< AphyloPruner > p(ptr);
  
  int ans = (int) sizeof(*p);
  
  Rcpp::Rcout << ans << std::endl;
  return ans;
}

// [[Rcpp::export(name = ".LogLike_pruner", rng = false)]]
List LogLike_pruner(
    SEXP tree_ptr,
    const std::vector< double > & mu_d,
    const std::vector< double > & mu_s,
    const std::vector< double > & psi,
    const std::vector< double > & eta,
    const double & Pi,
    bool verb = true,
    bool check_dims = false
) {
  Rcpp::XPtr< AphyloPruner > p(tree_ptr);
  
  // Setting the parameters
  p->args->set_mu_d(mu_d);
  p->args->set_mu_s(mu_s);
  p->args->set_psi(psi);
  
  p->args->set_eta(eta);
  
  // In the case of Pi, if it is negative, then it means that we are using
  // the stationary value of the transition probabilities.
  if (Pi < 0.0) {
    p->args->set_pi(
        (1 - p->args->prop_type_d)* mu_s[0]/(mu_s[0] + mu_s[1]) +
          p->args->prop_type_d * mu_d[0]/(mu_d[0] + mu_d[1])
    );
  } else 
    p->args->set_pi(Pi);
  
  // Calculating likelihood using Felsestein's algorithm.
  p->prune_postorder();
  
  if (verb) {
    NumericMatrix Pr(p->args->n, p->args->nstates);
    for (unsigned int i = 0u; i < p->args->n; ++i)
      for (unsigned int j = 0u; j < p->args->nstates; ++j)
        Pr(i, j) = p->args->Pr[i][j];
    
    return List::create(
      _["Pr"] = List::create(Pr),
      _["ll"] = wrap(p->args->ll)
      );
  } else
    return List::create(_["ll"] = wrap(p->args->ll));
}

// [[Rcpp::export(rng = false)]]
std::vector< std::vector< unsigned int > > Tree_get_offspring(const SEXP & tree_ptr) {
  
  Rcpp::XPtr< AphyloPruner > p(tree_ptr);
  return p->get_offspring();
  
}

// [[Rcpp::export(rng = false)]]
std::vector< std::vector< unsigned int > > Tree_get_parents(const SEXP & tree_ptr) {
  
  Rcpp::XPtr< AphyloPruner > p(tree_ptr);
  return p->get_parents();
  
}


// [[Rcpp::export(name=".Nnode_aphylo_pruner", rng = false)]]
unsigned int Tree_Nnode(const SEXP & tree_ptr, bool internal_only = true) {
  
  Rcpp::XPtr< AphyloPruner > p(tree_ptr);
  
  unsigned int count = p->n_nodes();
  
  if (internal_only)
    count -= p->n_tips();
  
  return count;
}

//' @rdname new_aphylo_pruner
//' @param ptr An object of class `aphylo_pruner`.
//' @return `dist2root`: An integer vector with the number of steps from each
//' node (internal or not) to the root node.
//' @export
// [[Rcpp::export(name = "dist2root", rng = false)]]
std::vector< unsigned int > Tree_get_dist_tip2root(const SEXP & ptr) {
  
  if (!Rf_inherits(ptr, "aphylo_pruner"))
    stop("-ptr- must be an object of class 'aphylo_pruner'.");
  
  Rcpp::XPtr< AphyloPruner > p(ptr);
  pruner::v_uint ans = p->get_dist_tip2root(), ans_sorted;
  pruner::v_uint tip = p->get_tips();
  
  // Right sorting
  ans_sorted.resize(ans.size());
  for (unsigned int i = 0u; i < tip.size(); ++i)
    ans_sorted[tip[i]] = ans[i];
  
  return ans_sorted;
  
}

//' @rdname new_aphylo_pruner
//' @return `get_postorder`: An integer vector with the postorder sequence
//' for pruning the tree (indexed from 0).
//' @export
// [[Rcpp::export(name = "get_postorder", rng = false)]]
std::vector< unsigned int > Tree_get_postorder(const SEXP & ptr) {
  
  if (!Rf_inherits(ptr, "aphylo_pruner"))
    stop("-ptr- must be an object of class 'aphylo_pruner'.");
  
  Rcpp::XPtr< AphyloPruner > p(ptr);
  
  return p->get_postorder();
  
}

//' @export
// [[Rcpp::export(name="Ntip.aphylo_pruner", rng = false)]]
unsigned int Tree_Ntip(const SEXP & phy) {
  
  Rcpp::XPtr< AphyloPruner > p(phy);
  
  return p->n_tips();
}


//' @export
// [[Rcpp::export(name="Nannotated.aphylo_pruner", rng = false)]]
unsigned int Tree_Nannotated(const SEXP & phy) {
  
  Rcpp::XPtr< AphyloPruner > p(phy);
  
  return p->args->nannotated;
}

//' @export
// [[Rcpp::export(name="Nann.aphylo_pruner", rng = false)]]
unsigned int Tree_Nann(const SEXP & phy) {
  
  Rcpp::XPtr< AphyloPruner > p(phy);
  
  return p->args->nfuns;
}


// [[Rcpp::export]]
unsigned int Tree_set_ann(const SEXP & phy, unsigned int i, unsigned int j, unsigned int val) {
  
  Rcpp::XPtr< AphyloPruner > p(phy);
  
  p->args->set_ann(i, j, val);
  
  return 0u;
  
}

// [[Rcpp::export]]
std::vector< std::vector< unsigned int > > Tree_get_ann(const SEXP & phy) {
  
  Rcpp::XPtr< AphyloPruner > p(phy);
  return p->args->A;
  
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
  types    = A
)

tree_ptr <- new_aphylo_pruner(dat)

aphylo_ll <- aphylo::LogLike

microbenchmark::microbenchmark(
  new = aphylo::LogLike(tree_ptr, mu_d = mu, mu_s=mu, psi = psi, eta = eta, Pi = Pi),
  old = aphylo::LogLike(dat, psi = psi, mu_d = mu, mu_s=mu, Pi = Pi, eta = eta)$ll,
  times = 100, unit = "relative"
)


*/
