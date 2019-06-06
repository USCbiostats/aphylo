#' Computes Log-likelihood
#' 
#' This function computes the log-likelihood of the chosen parameters given
#' a particular dataset. The arguments \code{annotations}, and \code{offspring}
#' should be as those returned by \code{\link{new_aphylo}}.
#' For complete Maximum Likelihood Estimation see [aphylo_estimates][aphylo_estimates-class].
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
#' \item{\code{mu}: A vector of length 2 with \eqn{\mu_0}{mu[0]} and
#' \eqn{\mu_1}{mu[1]} which are the gain and loss probabilities respectively.}
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
  mu,
  eta,
  Pi, 
  verb_ans    = TRUE,
  check_dims  = TRUE
) UseMethod("LogLike")

#' @export
#' @rdname LogLike
LogLike.aphylo <- function(
  tree,
  psi,
  mu,
  eta,
  Pi, 
  verb_ans    = TRUE,
  check_dims  = TRUE
) {
  
  ans <- .LogLike(
    annotations = list(with(tree, rbind(tip.annotation, node.annotation))),
    offspring   = list(tree$offspring),
    pseq        = list(tree$reduced_pseq),
    psi         = psi,
    mu          = mu,
    eta         = eta,
    Pi          = Pi,
    verb_ans    = verb_ans,
    check_dims  = check_dims
  )
  
  if (verb_ans)
    dim(ans$Pr[[1]]) <- c(
      ape::Nnode(tree, internal.only = FALSE),
      length(ans$Pr[[1]]) %/% ape::Nnode(tree, internal.only = FALSE)
    )
  
  ans
  
}

#' @export
#' @rdname LogLike
LogLike.multiAphylo <- function(
  tree,
  psi,
  mu,
  eta,
  Pi, 
  verb_ans    = TRUE,
  check_dims  = TRUE
) {
 
  annotations. <- lapply(tree, function(t.) with(t., rbind(tip.annotation, node.annotation)))
  offspring.   <- lapply(tree, "[[", "offspring")
  pseq.        <- lapply(tree, "[[", "reduced_pseq")
   
  ans <- .LogLike(
    annotations = annotations.,
    offspring   = offspring.,
    pseq        = pseq.,
    psi         = psi,
    mu          = mu,
    eta         = eta,
    Pi          = Pi,
    verb_ans    = verb_ans,
    check_dims  = check_dims
  )
  
  if (verb_ans)
    for (i in seq_along(ans$Pr))
    dim(ans$Pr[[i]]) <- c(
      ape::Nnode(tree[[i]], internal.only = FALSE),
      length(ans$Pr[[i]]) %/% ape::Nnode(tree[[i]], internal.only = FALSE)
    )
  
  ans
  
}