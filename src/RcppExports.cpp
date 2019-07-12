// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// auc
List auc(NumericVector pred, IntegerVector labels, int nc, bool nine_na);
RcppExport SEXP _aphylo_auc(SEXP predSEXP, SEXP labelsSEXP, SEXP ncSEXP, SEXP nine_naSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pred(predSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< bool >::type nine_na(nine_naSEXP);
    rcpp_result_gen = Rcpp::wrap(auc(pred, labels, nc, nine_na));
    return rcpp_result_gen;
END_RCPP
}
// new_aphylo_pruner
SEXP new_aphylo_pruner(const std::vector< std::vector< unsigned int > >& edgelist, const std::vector< std::vector< unsigned int > >& A, const std::vector< unsigned int >& Ntype);
RcppExport SEXP _aphylo_new_aphylo_pruner(SEXP edgelistSEXP, SEXP ASEXP, SEXP NtypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const std::vector< std::vector< unsigned int > >& >::type edgelist(edgelistSEXP);
    Rcpp::traits::input_parameter< const std::vector< std::vector< unsigned int > >& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const std::vector< unsigned int >& >::type Ntype(NtypeSEXP);
    rcpp_result_gen = Rcpp::wrap(new_aphylo_pruner(edgelist, A, Ntype));
    return rcpp_result_gen;
END_RCPP
}
// LogLike_pruner
List LogLike_pruner(SEXP tree_ptr, const std::vector< double >& mu, const std::vector< double >& psi, const std::vector< double >& eta, const double& Pi, bool verb, bool check_dims);
RcppExport SEXP _aphylo_LogLike_pruner(SEXP tree_ptrSEXP, SEXP muSEXP, SEXP psiSEXP, SEXP etaSEXP, SEXP PiSEXP, SEXP verbSEXP, SEXP check_dimsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type tree_ptr(tree_ptrSEXP);
    Rcpp::traits::input_parameter< const std::vector< double >& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const std::vector< double >& >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const std::vector< double >& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const double& >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< bool >::type verb(verbSEXP);
    Rcpp::traits::input_parameter< bool >::type check_dims(check_dimsSEXP);
    rcpp_result_gen = Rcpp::wrap(LogLike_pruner(tree_ptr, mu, psi, eta, Pi, verb, check_dims));
    return rcpp_result_gen;
END_RCPP
}
// Tree_get_offspring
std::vector< std::vector< unsigned int > > Tree_get_offspring(const SEXP& tree_ptr);
RcppExport SEXP _aphylo_Tree_get_offspring(SEXP tree_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type tree_ptr(tree_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(Tree_get_offspring(tree_ptr));
    return rcpp_result_gen;
END_RCPP
}
// Tree_get_parents
std::vector< std::vector< unsigned int > > Tree_get_parents(const SEXP& tree_ptr);
RcppExport SEXP _aphylo_Tree_get_parents(SEXP tree_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type tree_ptr(tree_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(Tree_get_parents(tree_ptr));
    return rcpp_result_gen;
END_RCPP
}
// Tree_Nnode
unsigned int Tree_Nnode(const SEXP& tree_ptr, bool internal_only);
RcppExport SEXP _aphylo_Tree_Nnode(SEXP tree_ptrSEXP, SEXP internal_onlySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type tree_ptr(tree_ptrSEXP);
    Rcpp::traits::input_parameter< bool >::type internal_only(internal_onlySEXP);
    rcpp_result_gen = Rcpp::wrap(Tree_Nnode(tree_ptr, internal_only));
    return rcpp_result_gen;
END_RCPP
}
// Tree_Ntip
unsigned int Tree_Ntip(const SEXP& phy);
RcppExport SEXP _aphylo_Tree_Ntip(SEXP phySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type phy(phySEXP);
    rcpp_result_gen = Rcpp::wrap(Tree_Ntip(phy));
    return rcpp_result_gen;
END_RCPP
}
// Tree_Nann
unsigned int Tree_Nann(const SEXP& phy);
RcppExport SEXP _aphylo_Tree_Nann(SEXP phySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type phy(phySEXP);
    rcpp_result_gen = Rcpp::wrap(Tree_Nann(phy));
    return rcpp_result_gen;
END_RCPP
}
// approx_geodesic
arma::imat approx_geodesic(const arma::umat& edges, unsigned int nsteps, bool undirected, bool warn);
RcppExport SEXP _aphylo_approx_geodesic(SEXP edgesSEXP, SEXP nstepsSEXP, SEXP undirectedSEXP, SEXP warnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::umat& >::type edges(edgesSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nsteps(nstepsSEXP);
    Rcpp::traits::input_parameter< bool >::type undirected(undirectedSEXP);
    Rcpp::traits::input_parameter< bool >::type warn(warnSEXP);
    rcpp_result_gen = Rcpp::wrap(approx_geodesic(edges, nsteps, undirected, warn));
    return rcpp_result_gen;
END_RCPP
}
// states
arma::imat states(int P);
RcppExport SEXP _aphylo_states(SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(states(P));
    return rcpp_result_gen;
END_RCPP
}
// prob_mat
arma::mat prob_mat(const arma::vec& pr);
RcppExport SEXP _aphylo_prob_mat(SEXP prSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type pr(prSEXP);
    rcpp_result_gen = Rcpp::wrap(prob_mat(pr));
    return rcpp_result_gen;
END_RCPP
}
// reduce_pseq
IntegerVector reduce_pseq(const arma::ivec& pseq, const arma::mat& A, const List& offspring);
RcppExport SEXP _aphylo_reduce_pseq(SEXP pseqSEXP, SEXP ASEXP, SEXP offspringSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::ivec& >::type pseq(pseqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const List& >::type offspring(offspringSEXP);
    rcpp_result_gen = Rcpp::wrap(reduce_pseq(pseq, A, offspring));
    return rcpp_result_gen;
END_RCPP
}
// root_node_prob
arma::vec root_node_prob(double Pi, const arma::imat& S);
RcppExport SEXP _aphylo_root_node_prob(SEXP PiSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(root_node_prob(Pi, S));
    return rcpp_result_gen;
END_RCPP
}
// probabilities
arma::mat probabilities(const arma::imat& annotations, const arma::ivec& pseq, const arma::vec& psi, const arma::vec& mu, const arma::vec& eta, const arma::imat& S, const ListOf<IntegerVector>& offspring);
RcppExport SEXP _aphylo_probabilities(SEXP annotationsSEXP, SEXP pseqSEXP, SEXP psiSEXP, SEXP muSEXP, SEXP etaSEXP, SEXP SSEXP, SEXP offspringSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::imat& >::type annotations(annotationsSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pseq(pseqSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const ListOf<IntegerVector>& >::type offspring(offspringSEXP);
    rcpp_result_gen = Rcpp::wrap(probabilities(annotations, pseq, psi, mu, eta, S, offspring));
    return rcpp_result_gen;
END_RCPP
}
// LogLike
List LogLike(const std::vector< arma::imat >& annotations, const std::vector< ListOf<IntegerVector> >& offspring, const std::vector< arma::ivec >& pseq, const arma::vec& psi, const arma::vec& mu, const arma::vec& eta, double Pi, bool verb_ans, bool check_dims);
RcppExport SEXP _aphylo_LogLike(SEXP annotationsSEXP, SEXP offspringSEXP, SEXP pseqSEXP, SEXP psiSEXP, SEXP muSEXP, SEXP etaSEXP, SEXP PiSEXP, SEXP verb_ansSEXP, SEXP check_dimsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const std::vector< arma::imat >& >::type annotations(annotationsSEXP);
    Rcpp::traits::input_parameter< const std::vector< ListOf<IntegerVector> >& >::type offspring(offspringSEXP);
    Rcpp::traits::input_parameter< const std::vector< arma::ivec >& >::type pseq(pseqSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< bool >::type verb_ans(verb_ansSEXP);
    Rcpp::traits::input_parameter< bool >::type check_dims(check_dimsSEXP);
    rcpp_result_gen = Rcpp::wrap(LogLike(annotations, offspring, pseq, psi, mu, eta, Pi, verb_ans, check_dims));
    return rcpp_result_gen;
END_RCPP
}
// posterior_prob
List posterior_prob(const arma::mat& Pr_postorder, const arma::vec& mu, const double& Pi, const arma::ivec& pseq, const List& offspring);
RcppExport SEXP _aphylo_posterior_prob(SEXP Pr_postorderSEXP, SEXP muSEXP, SEXP PiSEXP, SEXP pseqSEXP, SEXP offspringSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Pr_postorder(Pr_postorderSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const double& >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pseq(pseqSEXP);
    Rcpp::traits::input_parameter< const List& >::type offspring(offspringSEXP);
    rcpp_result_gen = Rcpp::wrap(posterior_prob(Pr_postorder, mu, Pi, pseq, offspring));
    return rcpp_result_gen;
END_RCPP
}
// prediction_score_rand
double prediction_score_rand(const arma::mat& A, const arma::mat& W, double alpha);
RcppExport SEXP _aphylo_prediction_score_rand(SEXP ASEXP, SEXP WSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(prediction_score_rand(A, W, alpha));
    return rcpp_result_gen;
END_RCPP
}
// sim_fun_on_tree
IntegerMatrix sim_fun_on_tree(const List& offspring, const arma::ivec& pseq, const arma::vec& psi, const arma::vec& mu, const arma::vec& eta, const arma::vec& Pi, int P);
RcppExport SEXP _aphylo_sim_fun_on_tree(SEXP offspringSEXP, SEXP pseqSEXP, SEXP psiSEXP, SEXP muSEXP, SEXP etaSEXP, SEXP PiSEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type offspring(offspringSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type pseq(pseqSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< int >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_fun_on_tree(offspring, pseq, psi, mu, eta, Pi, P));
    return rcpp_result_gen;
END_RCPP
}
// sim_tree
List sim_tree(int n, Function f, bool branches);
RcppExport SEXP _aphylo_sim_tree(SEXP nSEXP, SEXP fSEXP, SEXP branchesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< Function >::type f(fSEXP);
    Rcpp::traits::input_parameter< bool >::type branches(branchesSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_tree(n, f, branches));
    return rcpp_result_gen;
END_RCPP
}
// fast_table
IntegerMatrix fast_table(const arma::ivec& x);
RcppExport SEXP _aphylo_fast_table(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::ivec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_table(x));
    return rcpp_result_gen;
END_RCPP
}
// fast_table_using_labels
arma::uvec fast_table_using_labels(const arma::ivec& x, const arma::ivec& ids);
RcppExport SEXP _aphylo_fast_table_using_labels(SEXP xSEXP, SEXP idsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::ivec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ids(idsSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_table_using_labels(x, ids));
    return rcpp_result_gen;
END_RCPP
}
// list_offspring
ListOf<IntegerVector> list_offspring(IntegerMatrix E, int n);
RcppExport SEXP _aphylo_list_offspring(SEXP ESEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type E(ESEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(list_offspring(E, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_aphylo_auc", (DL_FUNC) &_aphylo_auc, 4},
    {"_aphylo_new_aphylo_pruner", (DL_FUNC) &_aphylo_new_aphylo_pruner, 3},
    {"_aphylo_LogLike_pruner", (DL_FUNC) &_aphylo_LogLike_pruner, 7},
    {"_aphylo_Tree_get_offspring", (DL_FUNC) &_aphylo_Tree_get_offspring, 1},
    {"_aphylo_Tree_get_parents", (DL_FUNC) &_aphylo_Tree_get_parents, 1},
    {"_aphylo_Tree_Nnode", (DL_FUNC) &_aphylo_Tree_Nnode, 2},
    {"_aphylo_Tree_Ntip", (DL_FUNC) &_aphylo_Tree_Ntip, 1},
    {"_aphylo_Tree_Nann", (DL_FUNC) &_aphylo_Tree_Nann, 1},
    {"_aphylo_approx_geodesic", (DL_FUNC) &_aphylo_approx_geodesic, 4},
    {"_aphylo_states", (DL_FUNC) &_aphylo_states, 1},
    {"_aphylo_prob_mat", (DL_FUNC) &_aphylo_prob_mat, 1},
    {"_aphylo_reduce_pseq", (DL_FUNC) &_aphylo_reduce_pseq, 3},
    {"_aphylo_root_node_prob", (DL_FUNC) &_aphylo_root_node_prob, 2},
    {"_aphylo_probabilities", (DL_FUNC) &_aphylo_probabilities, 7},
    {"_aphylo_LogLike", (DL_FUNC) &_aphylo_LogLike, 9},
    {"_aphylo_posterior_prob", (DL_FUNC) &_aphylo_posterior_prob, 5},
    {"_aphylo_prediction_score_rand", (DL_FUNC) &_aphylo_prediction_score_rand, 3},
    {"_aphylo_sim_fun_on_tree", (DL_FUNC) &_aphylo_sim_fun_on_tree, 7},
    {"_aphylo_sim_tree", (DL_FUNC) &_aphylo_sim_tree, 3},
    {"_aphylo_fast_table", (DL_FUNC) &_aphylo_fast_table, 1},
    {"_aphylo_fast_table_using_labels", (DL_FUNC) &_aphylo_fast_table_using_labels, 2},
    {"_aphylo_list_offspring", (DL_FUNC) &_aphylo_list_offspring, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_aphylo(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
