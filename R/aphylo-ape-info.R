#' APE wrappers
#' @param phy An object of class [aphylo] (see [ape::Nnode]).
#' @param ... Further arguments passed to the method.
#' @details 
#' Functions `Nedge`, `Nnode`, and `Ntip` have S3 methods in the `ape` package.
#' 
#' @return Integer with the number of edges, nodes, or tips accordignly.
#' @name ape-methods
#' @family information
NULL

#' @rdname ape-methods
#' @export
Nedge.aphylo <- function(phy) ape::Nedge(phy$tree)

#' @export
Nedge.aphylo_estimates <- function(phy) ape::Nedge(phy$dat)

#' @export
Nedge.multiAphylo <- function(phy) sapply(phy, ape::Nedge)

#' @export
Nedge.multiAphylo_pruner <- Nedge.multiAphylo

#' @export
Nnode.aphylo <- function(phy, ...) ape::Nnode(phy$tree, ...)

#'@export
Nnode.aphylo_estimates <- function(phy, ...) ape::Nnode(phy$dat, ...)

#' @export
Nnode.multiAphylo <- function(phy, ...) sapply(phy, ape::Nnode, ...)

#' @export
Nnode.aphylo_pruner <- function(phy, ...) {
  
  dots <- list(...)
  if ("internal.only" %in% names(dots))
    .Nnode_aphylo_pruner(phy, internal_only = dots$internal.only)
  else
    .Nnode_aphylo_pruner(phy, internal_only = TRUE)
  
}

#' @export
Nnode.multiAphylo_pruner <- Nnode.multiAphylo

#' @export
Ntip.aphylo <- function(phy) ape::Ntip(phy$tree)

#' @export
Ntip.aphylo_estimates <- function(phy) ape::Ntip(phy$dat)

#' @export
Ntip.multiAphylo <- function(phy) sapply(phy, ape::Ntip)

#' @export
Ntip.multiAphylo_pruner <- function(phy) {
  sapply(phy, ape::Ntip)
}

#' Information about `aphylo` and `multiAphylo` objects
#' 
#' Information about annotations, in particular, number of annotations (`Nann`),
#' number of annotated leaves (`Nannotated`), number of unnanotated leaves
#' (`Nunannotated`), and number of trees (`Ntrees`).
#' 
#' @param phy Either an object of class [aphylo], [multiAphylo], or
#' [aphylo_estimates].
#' 
#' @return If `phy` is of class `aphylo`, then a single scalar.
#' otherwise, if `phy` is of class `multiAphylo`
#' 
#' @name aphylo-info
#' @family information
NULL

#' @rdname aphylo-info
#' @export
Nann <- function(phy) UseMethod("Nann")

#' @export
Nann.aphylo <- function(phy) ncol(phy$tip.annotation)

#' @export
Nann.multiAphylo <- function(phy) sapply(phy, Nann)

#' @export
Nann.multiAphylo_pruner <- Nann.multiAphylo

#' @export
Nann.aphylo_estimates <- function(phy) {
  Nann(phy$dat)
}

#' @export
#' @rdname aphylo-info
Nannotated <- function(phy) UseMethod("Nannotated")

#' @export
Nannotated.aphylo <- function(phy) phy$Ntips.annotated

#' @export
Nannotated.multiAphylo <- function(phy) sapply(phy, Nannotated)

#' @export
Nannotated.multiAphylo_pruner <- Nannotated.multiAphylo


#' @export
Nannotated.aphylo_estimates <- function(phy) {
  Nannotated(phy$dat)
}

#' @export
#' @rdname aphylo-info
Ntrees <- function(phy) UseMethod("Ntrees")

#' @export
Ntrees.aphylo <- function(phy) 1L

#' @export
Ntrees.phylo <- function(phy) 1L

#' @export
Ntrees.multiPhylo <- function(phy) length(phy)

#' @export
Ntrees.multiAphylo <- function(phy) length(phy)

#' @export
Ntrees.multiAphylo_pruner <- Ntrees.multiAphylo

#' @export
Ntrees.aphylo_pruner <- Ntrees.aphylo

#' @export
Ntrees.aphylo_estimates <- function(phy) Ntrees(phy$dat)
