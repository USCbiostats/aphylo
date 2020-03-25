#' Likelihood of an observed annotated phylogenetic tree
#' 
#' This function computes the log-likelihood of the chosen parameters given
#' a particular dataset. The arguments `annotations`, and `offspring`
#' should be as those returned by [new_aphylo()].
#' For complete parameter estimation see [aphylo_estimates].
#' 
#' @template parameters
#' @param tree A phylogenetic tree of class [aphylo].
#' @templateVar .psi 1
#' @templateVar .mu 1
#' @templateVar .eta 1
#' @templateVar .Pi 1
#' @param verb_ans Logical scalar. When \code{FALSE} (default) the function
#' returns a list with a single scalar (the log-likelihood).
#' @param check_dims Logical scalar. When \code{TRUE} (default) the function
#' checks the dimmension of the passed parameters.
#' 
#' @details
#' The parameters to estimate are described as follows:
#' \enumerate{
#' \item{\code{psi}: A vector of length 2 with \eqn{\psi_0}{psi[0]} and
#' \eqn{\psi_1}{psi[1]}, which are the misclassification probabilities fo
#' \eqn{s_p=0}{s[p]=0} and \eqn{s_p=1}{s[p]=1}
#' respectively.}
#' \item{\code{mu_d}, \code{mu_s}: A vector of length 2 with \eqn{\mu_0}{mu[0]} and
#' \eqn{\mu_1}{mu[1]} which are the gain and loss probabilities respectively.
#' The subscript d denotes duplication nodes and s speciation node.}
#' \item{\code{eta}: A vector of length 2 with \eqn{\eta_0}{eta[0]} and
#' \eqn{\eta_1}{eta[1]} which are the annotation bias probabilities.}
#' \item{\code{Pi}: A numeric scalar which for which equals the probability
#' of the root node having the function.}
#' }
#' @return A list of class \code{phylo_LogLik} with the following elements:
#' \item{S}{An integer matrix of size \eqn{2^p\times p}{2^p * p} as returned
#' by \code{\link{states}}.}
#' \item{Pr}{A numeric matrix of size \eqn{G\times 2^p}{G * 2^p} with node/state
#' probabilities.}
#' \item{ll}{A numeric scalar with the log-likelihood value given the chosen
#' parameters.}
#' @export
LogLike <- function(
  tree,
  psi,
  mu_d,
  mu_s,
  eta,
  Pi, 
  verb_ans    = TRUE,
  check_dims  = TRUE
) UseMethod("LogLike")

#' @export
LogLike.aphylo_pruner <- function(
  tree,
  psi,
  mu_d,
  mu_s,
  eta,
  Pi, 
  verb_ans    = TRUE,
  check_dims  = TRUE
) {
  
  .LogLike_pruner(
    tree_ptr = tree,
    mu_d      = mu_d,
    mu_s      = mu_s,
    psi      = psi,
    eta      = eta,
    Pi       = Pi,
    verb     = verb_ans
  )
  
}

#' @export
LogLike.aphylo <- function(
  tree,
  psi,
  mu_d,
  mu_s,
  eta,
  Pi, 
  verb_ans    = TRUE,
  check_dims  = TRUE
) {
  
  tree_ptr <- new_aphylo_pruner(tree)
  .LogLike_pruner(
    tree_ptr = tree_ptr,
    mu_d      = mu_d,
    mu_s      = mu_s,
    psi      = psi,
    eta      = eta,
    Pi       = Pi,
    verb     = verb_ans
  )
  
}

#' @export
LogLike.multiAphylo <- function(
  tree,
  psi,
  mu_d,
  mu_s,
  eta,
  Pi, 
  verb_ans    = TRUE,
  check_dims  = TRUE
) {
 
  res <- list(ll = 0.0, Pr = NULL)
  for (i in seq_along(tree)) {
    
    tmp <- LogLike(
      tree[[i]],
      mu_d     = mu_d,
      mu_s     = mu_s,
      psi      = psi,
      eta      = eta,
      Pi       = Pi,
      verb_ans = verb_ans
      )
    
    res$ll <- res$ll + tmp$ll
    if (verb_ans)
      res$Pr <- c(res$Pr, tmp$Pr)
    
  }
  res
  
}

#' @export
LogLike.multiAphylo_pruner <- LogLike.multiAphylo

