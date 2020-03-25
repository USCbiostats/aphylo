#' @param params A vector of length 7 with initial parameters. In particular
#' `psi[1]`, `psi[2]`, `mu[1]`, `mu[2]`, `eta[1]`, `eta[2]` and `Pi`.
#' @param model A model as specified in [aphylo-model].
#' @param priors A function to be used as prior for the model (see [bprior]).
#' @param check_informative Logical scalar. When `TRUE` the algorithm
#' stops with an error when the annotations are uninformative (either 0s or 1s).
#' @param reduced_pseq Logical. When `TRUE` it will use a reduced peeling sequence
#' in which it drops unannotated leafs. If the model includes `eta` this is set
#' to `FALSE`.
#' 
#' @return An object of class [aphylo_estimates].
NULL
