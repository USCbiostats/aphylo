#' Pointer to `pruner`
#' 
#' Creates an external pointer to an object of class `aphylo_pruner`. This is mostly
#' used to compute the model's likelihood function faster by reusing underlying
#' C++ class objects to store probability matrices and data. This is intended
#' for internal use only
#' 
#' @param edgelist a List two integer vectors.
#' @param types An integer vector of types of size `N`.
#' @aliases aphylo_pruner
#' 
#' @details The underlying implementation of the pruning function is based on the
#' pruner C++ library that implements Felsenstein's tree pruning algorithm.
#' See ttps://github.com/USCbiostats/pruner.
#' 
#' @examples
#' set.seed(1)
#' x  <- raphylo(10) 
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
#' @name new_aphylo_pruner
NULL

#' @export
#' @rdname new_aphylo_pruner
#' @param ... Further arguments passed to the method
new_aphylo_pruner <- function(...) UseMethod("new_aphylo_pruner")

#' @export
#' @param annotation a list of length `N` (annotations).
#' @rdname new_aphylo_pruner
new_aphylo_pruner.matrix <- function(edgelist, annotation, types = NULL, ...) {
  
  n <- max(edgelist)
  P <- ncol(annotation)
  
  if (ncol(edgelist) != 2L)
    stop("`edgelist` must be a two-column matrix.", call. = FALSE)
  
  if (nrow(annotation) != n)
    stop("`nrow(annotation) != nrow(edge`", call. = FALSE)
  
  if (is.null(types))
    types <- rep(0L, n)
  else if (length(types) != n)
    stop("`nrow(annotation) != nrow(edge`", call. = FALSE)
    
  # Checking the complete cases
  nannotated <- annotation
  nannotated[nannotated == 9L] <- NA
  nannotated <- sum(complete.cases(nannotated))
  
  # Splitting the data
  annotation <- lapply(seq_len(n), function(i) as.vector(annotation[i, ]))
  
  new_aphylo_pruner.(
    edgelist = list(edgelist[, 1L] - 1L, edgelist[, 2L] - 1L),
    A        = annotation,
    types    = types
    )
  
}

#' @export
#' @rdname new_aphylo_pruner
#' @param x An object of class [aphylo].
new_aphylo_pruner.aphylo <- function(x, ...) {
  
  annotation <- with(x, rbind(tip.annotation, node.annotation))
  annotation <- lapply(seq_len(nrow(annotation)), function(i) annotation[i, ])
  
  new_aphylo_pruner.(
    edgelist   = list(x$tree$edge[, 1L] - 1L, x$tree$edge[, 2L] - 1L),
    A          = annotation,
    types      = x$types, 
    nannotated = x$Ntips.annotated
  )
  
}

#' @export
#' @rdname new_aphylo_pruner
new_aphylo_pruner.multiAphylo <- function(x, ...) {
  
  structure(
    lapply(x, new_aphylo_pruner, ...),
    class = "multiAphylo_pruner"
  )
  
}

