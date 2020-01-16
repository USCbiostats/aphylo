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

#' @rdname ape-methods
#' @export
Nedge.aphylo_estimates <- function(phy) ape::Nedge(phy$dat)

#' @rdname ape-methods
#' @export
Nedge.multiAphylo <- function(phy) sapply(phy, ape::Nedge)

#' @rdname ape-methods
#' @export
Nedge.multiAphylo_pruner <- Nedge.multiAphylo

#' @rdname ape-methods
#' @export
Nnode.aphylo <- function(phy, ...) ape::Nnode(phy$tree, ...)

#'@export
#'@rdname ape-methods
Nnode.aphylo_estimates <- function(phy, ...) ape::Nnode(phy$dat, ...)

#' @rdname ape-methods
#' @export
Nnode.multiAphylo <- function(phy, ...) sapply(phy, ape::Nnode, ...)

#' @rdname ape-methods
#' @export
Nnode.aphylo_pruner <- function(phy, ...) {
  
  dots <- list(...)
  if ("internal.only" %in% names(dots))
    .Nnode_aphylo_pruner(phy, internal_only = dots$internal.only)
  else
    .Nnode_aphylo_pruner(phy, internal_only = TRUE)
  
}

#' @rdname ape-methods
#' @export
Nnode.multiAphylo_pruner <- Nnode.multiAphylo

#' @rdname ape-methods
#' @export
Ntip.aphylo <- function(phy) ape::Ntip(phy$tree)

#' @export
#' @rdname ape-methods
Ntip.aphylo_estimates <- function(phy) ape::Ntip(phy$dat)

#' @export
#' @rdname ape-methods
Ntip.multiAphylo <- function(phy) sapply(phy, ape::Ntip)

#' @export
#' @rdname ape-methods
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

#' @rdname aphylo-info
#' @export
Nann.aphylo <- function(phy) ncol(phy$tip.annotation)

#' @rdname aphylo-info
#' @export
Nann.multiAphylo <- function(phy) sapply(phy, Nann)

#' @rdname aphylo-info
#' @export
Nann.multiAphylo_pruner <- Nann.multiAphylo

#' @export
#' @rdname aphylo-info
Nann.aphylo_estimates <- function(phy) {
  Nann(phy$dat)
}



#' @export
#' @rdname aphylo-info
Nannotated <- function(phy) UseMethod("Nannotated")

#' @export
#' @rdname aphylo-info
Nannotated.aphylo <- function(phy) phy$Ntips.annotated

#' @export
#' @rdname aphylo-info
Nannotated.multiAphylo <- function(phy) sapply(phy, Nannotated)

#' @export
#' @rdname aphylo-info
Nannotated.multiAphylo_pruner <- Nannotated.multiAphylo


#' @export
#' @rdname aphylo-info
Nannotated.aphylo_estimates <- function(phy) {
  Nannotated(phy$dat)
}

#' @export
#' @rdname aphylo-info
Ntrees <- function(phy) UseMethod("Ntrees")

#' @export
#' @rdname aphylo-info
Ntrees.aphylo <- function(phy) 1L

#' @export
#' @rdname aphylo-info
Ntrees.phylo <- function(phy) 1L

#' @export
#' @rdname aphylo-info
Ntrees.multiPhylo <- function(phy) length(phy)

#' @export
#' @rdname aphylo-info
Ntrees.multiAphylo <- function(phy) length(phy)

#' @export
#' @rdname aphylo-info
Ntrees.multiAphylo_pruner <- Ntrees.multiAphylo

#' @export
#' @rdname aphylo-info
Ntrees.aphylo_pruner <- Ntrees.aphylo

#' @export
#' @rdname aphylo-info
Ntrees.aphylo_estimates <- function(phy) Ntrees(phy$dat)
