// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// new_aphylo_pruner_cpp
SEXP new_aphylo_pruner_cpp(const std::vector< std::vector< unsigned int > >& edgelist, const std::vector< std::vector< unsigned int > >& A, const std::vector< unsigned int >& types, unsigned int nannotated);
RcppExport SEXP _aphylo_new_aphylo_pruner_cpp(SEXP edgelistSEXP, SEXP ASEXP, SEXP typesSEXP, SEXP nannotatedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const std::vector< std::vector< unsigned int > >& >::type edgelist(edgelistSEXP);
    Rcpp::traits::input_parameter< const std::vector< std::vector< unsigned int > >& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const std::vector< unsigned int >& >::type types(typesSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nannotated(nannotatedSEXP);
    rcpp_result_gen = Rcpp::wrap(new_aphylo_pruner_cpp(edgelist, A, types, nannotated));
    return rcpp_result_gen;
END_RCPP
}
// sizeof_pruner
int sizeof_pruner(SEXP ptr);
RcppExport SEXP _aphylo_sizeof_pruner(SEXP ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type ptr(ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(sizeof_pruner(ptr));
    return rcpp_result_gen;
END_RCPP
}
// LogLike_pruner
List LogLike_pruner(SEXP tree_ptr, const std::vector< double >& mu_d, const std::vector< double >& mu_s, const std::vector< double >& psi, const std::vector< double >& eta, const double& Pi, bool verb, bool check_dims);
RcppExport SEXP _aphylo_LogLike_pruner(SEXP tree_ptrSEXP, SEXP mu_dSEXP, SEXP mu_sSEXP, SEXP psiSEXP, SEXP etaSEXP, SEXP PiSEXP, SEXP verbSEXP, SEXP check_dimsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type tree_ptr(tree_ptrSEXP);
    Rcpp::traits::input_parameter< const std::vector< double >& >::type mu_d(mu_dSEXP);
    Rcpp::traits::input_parameter< const std::vector< double >& >::type mu_s(mu_sSEXP);
    Rcpp::traits::input_parameter< const std::vector< double >& >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const std::vector< double >& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const double& >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< bool >::type verb(verbSEXP);
    Rcpp::traits::input_parameter< bool >::type check_dims(check_dimsSEXP);
    rcpp_result_gen = Rcpp::wrap(LogLike_pruner(tree_ptr, mu_d, mu_s, psi, eta, Pi, verb, check_dims));
    return rcpp_result_gen;
END_RCPP
}
// Tree_get_offspring
std::vector< std::vector< unsigned int > > Tree_get_offspring(const SEXP& tree_ptr);
RcppExport SEXP _aphylo_Tree_get_offspring(SEXP tree_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
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
    Rcpp::traits::input_parameter< const SEXP& >::type tree_ptr(tree_ptrSEXP);
    Rcpp::traits::input_parameter< bool >::type internal_only(internal_onlySEXP);
    rcpp_result_gen = Rcpp::wrap(Tree_Nnode(tree_ptr, internal_only));
    return rcpp_result_gen;
END_RCPP
}
// Tree_get_dist_tip2root
std::vector< unsigned int > Tree_get_dist_tip2root(const SEXP& ptr);
RcppExport SEXP _aphylo_Tree_get_dist_tip2root(SEXP ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type ptr(ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(Tree_get_dist_tip2root(ptr));
    return rcpp_result_gen;
END_RCPP
}
// Tree_get_postorder
std::vector< unsigned int > Tree_get_postorder(const SEXP& ptr);
RcppExport SEXP _aphylo_Tree_get_postorder(SEXP ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type ptr(ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(Tree_get_postorder(ptr));
    return rcpp_result_gen;
END_RCPP
}
// Tree_Ntip
unsigned int Tree_Ntip(const SEXP& phy);
RcppExport SEXP _aphylo_Tree_Ntip(SEXP phySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type phy(phySEXP);
    rcpp_result_gen = Rcpp::wrap(Tree_Ntip(phy));
    return rcpp_result_gen;
END_RCPP
}
// Tree_Nannotated
unsigned int Tree_Nannotated(const SEXP& phy);
RcppExport SEXP _aphylo_Tree_Nannotated(SEXP phySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type phy(phySEXP);
    rcpp_result_gen = Rcpp::wrap(Tree_Nannotated(phy));
    return rcpp_result_gen;
END_RCPP
}
// Tree_Nann
unsigned int Tree_Nann(const SEXP& phy);
RcppExport SEXP _aphylo_Tree_Nann(SEXP phySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type phy(phySEXP);
    rcpp_result_gen = Rcpp::wrap(Tree_Nann(phy));
    return rcpp_result_gen;
END_RCPP
}
// Tree_set_ann
unsigned int Tree_set_ann(const SEXP& phy, unsigned int i, unsigned int j, unsigned int val);
RcppExport SEXP _aphylo_Tree_set_ann(SEXP phySEXP, SEXP iSEXP, SEXP jSEXP, SEXP valSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type phy(phySEXP);
    Rcpp::traits::input_parameter< unsigned int >::type i(iSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type j(jSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type val(valSEXP);
    rcpp_result_gen = Rcpp::wrap(Tree_set_ann(phy, i, j, val));
    return rcpp_result_gen;
END_RCPP
}
// Tree_get_ann
std::vector< std::vector< unsigned int > > Tree_get_ann(const SEXP& phy);
RcppExport SEXP _aphylo_Tree_get_ann(SEXP phySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type phy(phySEXP);
    rcpp_result_gen = Rcpp::wrap(Tree_get_ann(phy));
    return rcpp_result_gen;
END_RCPP
}
// auc
List auc(const NumericVector& pred, const IntegerVector& labels, int nc, bool nine_na);
RcppExport SEXP _aphylo_auc(SEXP predSEXP, SEXP labelsSEXP, SEXP ncSEXP, SEXP nine_naSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type pred(predSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< bool >::type nine_na(nine_naSEXP);
    rcpp_result_gen = Rcpp::wrap(auc(pred, labels, nc, nine_na));
    return rcpp_result_gen;
END_RCPP
}
// states
IntegerMatrix states(int P);
RcppExport SEXP _aphylo_states(SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< int >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(states(P));
    return rcpp_result_gen;
END_RCPP
}
// prob_mat
NumericMatrix prob_mat(const NumericVector& pr);
RcppExport SEXP _aphylo_prob_mat(SEXP prSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type pr(prSEXP);
    rcpp_result_gen = Rcpp::wrap(prob_mat(pr));
    return rcpp_result_gen;
END_RCPP
}
// root_node_prob
NumericVector root_node_prob(double Pi, const IntegerMatrix& S);
RcppExport SEXP _aphylo_root_node_prob(SEXP PiSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< double >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(root_node_prob(Pi, S));
    return rcpp_result_gen;
END_RCPP
}
// reduce_pseq
IntegerVector reduce_pseq(const IntegerVector& pseq, const NumericMatrix& A, const List& offspring);
RcppExport SEXP _aphylo_reduce_pseq(SEXP pseqSEXP, SEXP ASEXP, SEXP offspringSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type pseq(pseqSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const List& >::type offspring(offspringSEXP);
    rcpp_result_gen = Rcpp::wrap(reduce_pseq(pseq, A, offspring));
    return rcpp_result_gen;
END_RCPP
}
// posterior_prob
List posterior_prob(const NumericMatrix& Pr_postorder, const std::vector< unsigned int >& types, const NumericVector& mu_d, const NumericVector& mu_s, const double& Pi, const IntegerVector& pseq, const List& offspring);
RcppExport SEXP _aphylo_posterior_prob(SEXP Pr_postorderSEXP, SEXP typesSEXP, SEXP mu_dSEXP, SEXP mu_sSEXP, SEXP PiSEXP, SEXP pseqSEXP, SEXP offspringSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Pr_postorder(Pr_postorderSEXP);
    Rcpp::traits::input_parameter< const std::vector< unsigned int >& >::type types(typesSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type mu_d(mu_dSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type mu_s(mu_sSEXP);
    Rcpp::traits::input_parameter< const double& >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type pseq(pseqSEXP);
    Rcpp::traits::input_parameter< const List& >::type offspring(offspringSEXP);
    rcpp_result_gen = Rcpp::wrap(posterior_prob(Pr_postorder, types, mu_d, mu_s, Pi, pseq, offspring));
    return rcpp_result_gen;
END_RCPP
}
// sim_fun_on_tree
IntegerMatrix sim_fun_on_tree(const List& offspring, const IntegerVector& types, const IntegerVector& pseq, const NumericVector& psi, const NumericVector& mu_d, const NumericVector& mu_s, const NumericVector& eta, const NumericVector& Pi, int P);
RcppExport SEXP _aphylo_sim_fun_on_tree(SEXP offspringSEXP, SEXP typesSEXP, SEXP pseqSEXP, SEXP psiSEXP, SEXP mu_dSEXP, SEXP mu_sSEXP, SEXP etaSEXP, SEXP PiSEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type offspring(offspringSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type types(typesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type pseq(pseqSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type mu_d(mu_dSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type mu_s(mu_sSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< int >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_fun_on_tree(offspring, types, pseq, psi, mu_d, mu_s, eta, Pi, P));
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

static const R_CallMethodDef CallEntries[] = {
    {"_aphylo_new_aphylo_pruner_cpp", (DL_FUNC) &_aphylo_new_aphylo_pruner_cpp, 4},
    {"_aphylo_sizeof_pruner", (DL_FUNC) &_aphylo_sizeof_pruner, 1},
    {"_aphylo_LogLike_pruner", (DL_FUNC) &_aphylo_LogLike_pruner, 8},
    {"_aphylo_Tree_get_offspring", (DL_FUNC) &_aphylo_Tree_get_offspring, 1},
    {"_aphylo_Tree_get_parents", (DL_FUNC) &_aphylo_Tree_get_parents, 1},
    {"_aphylo_Tree_Nnode", (DL_FUNC) &_aphylo_Tree_Nnode, 2},
    {"_aphylo_Tree_get_dist_tip2root", (DL_FUNC) &_aphylo_Tree_get_dist_tip2root, 1},
    {"_aphylo_Tree_get_postorder", (DL_FUNC) &_aphylo_Tree_get_postorder, 1},
    {"_aphylo_Tree_Ntip", (DL_FUNC) &_aphylo_Tree_Ntip, 1},
    {"_aphylo_Tree_Nannotated", (DL_FUNC) &_aphylo_Tree_Nannotated, 1},
    {"_aphylo_Tree_Nann", (DL_FUNC) &_aphylo_Tree_Nann, 1},
    {"_aphylo_Tree_set_ann", (DL_FUNC) &_aphylo_Tree_set_ann, 4},
    {"_aphylo_Tree_get_ann", (DL_FUNC) &_aphylo_Tree_get_ann, 1},
    {"_aphylo_auc", (DL_FUNC) &_aphylo_auc, 4},
    {"_aphylo_states", (DL_FUNC) &_aphylo_states, 1},
    {"_aphylo_prob_mat", (DL_FUNC) &_aphylo_prob_mat, 1},
    {"_aphylo_root_node_prob", (DL_FUNC) &_aphylo_root_node_prob, 2},
    {"_aphylo_reduce_pseq", (DL_FUNC) &_aphylo_reduce_pseq, 3},
    {"_aphylo_posterior_prob", (DL_FUNC) &_aphylo_posterior_prob, 7},
    {"_aphylo_sim_fun_on_tree", (DL_FUNC) &_aphylo_sim_fun_on_tree, 9},
    {"_aphylo_sim_tree", (DL_FUNC) &_aphylo_sim_tree, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_aphylo(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
