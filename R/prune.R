#' Pointer to `pruner`
#' 
#' Creates an external pointer to an object of class `aphylo_pruner`. This is mostly
#' used to compute the model's likelihood function faster by reusing underlying
#' C++ class objects to store probability matrices and data. This is intended
#' for internal use only.
#' @param x An object of class [aphylo] or [multiAphylo].
#' @aliases aphylo_pruner
#' @return The function `new_aphylo_pruner` returns an object of class
#' `aphylo_pruner` or `multiAphylo_pruner`, depending on the class of `x`.
#' @details The underlying implementation of the pruning function is based on the
#' pruner C++ library that implements Felsenstein's tree pruning algorithm.
#' See \url{https://github.com/USCbiostats/pruner}.
#' 
#' @examples
#' set.seed(1)
#' x  <- raphylo(20) 
#' pruner <- new_aphylo_pruner(x)
#' 
#' # Computing loglike
#' LogLike(
#'   pruner,
#'   psi  = c(.10, .20),
#'   mu_d = c(.90, .80),
#'   mu_s = c(.10, .05),
#'   Pi   = .05,
#'   eta  = c(.90, .80)
#'   )
#'   
#' dist2root(pruner)
#' get_postorder(pruner)
#' @name new_aphylo_pruner
NULL

#' @export
#' @rdname new_aphylo_pruner
#' @param ... Further arguments passed to the method
new_aphylo_pruner <- function(x, ...) UseMethod("new_aphylo_pruner")

# new_aphylo_pruner.aphylo_pruner <- function(x, ...) x
# new_aphylo_pruner.multiAphylo_pruner <- function(x, ...) x

#' @export
new_aphylo_pruner.aphylo <- function(x, ...) {
  
  annotation <- with(x, rbind(tip.annotation, node.annotation))
  annotation <- lapply(seq_len(nrow(annotation)), function(i) annotation[i, ])
  
  new_aphylo_pruner_cpp(
    edgelist   = list(x$tree$edge[, 1L] - 1L, x$tree$edge[, 2L] - 1L),
    A          = annotation,
    types      = with(x, c(tip.type, node.type)), 
    nannotated = x$Ntips.annotated
  )
  
}

#' @export
new_aphylo_pruner.multiAphylo <- function(x, ...) {
  
  structure(
    lapply(x, new_aphylo_pruner, ...),
    class = "multiAphylo_pruner"
  )
  
}



