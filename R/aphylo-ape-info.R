#' Available methods from the APE package
#' 
#' The generics [ape::Nedge()], [ape::Nnode()], and [ape::Ntip()] can be used
#' directly on objects of class [aphylo], [aphylo_estimates], [multiAphylo]
#' 
#' @return Integer with the number of edges, nodes, or tips accordignly.
#' @name ape-methods
#' @family information
#' @importFrom ape as.phylo plot.phylo Nedge Ntip Nnode
#' @examples 
#' set.seed(12312)
#' atree <- raphylo(50, P = 2)
#' Nnode(atree)
#' Ntip(atree)
#' Nedge(atree)
#' 
#' multitree <- rmultiAphylo(10, 50, P = 2)
#' Nnode(multitree)
#' Ntip(multitree)
#' Nedge(multitree)
NULL

#' @method Nedge aphylo
#' @export
Nedge.aphylo <- function(phy) ape::Nedge(phy$tree)

#' @method Nedge aphylo_estimates
#' @export
Nedge.aphylo_estimates <- function(phy) ape::Nedge(phy$dat)

#' @method Nedge multiAphylo
#' @export
Nedge.multiAphylo <- function(phy) sapply(phy, ape::Nedge)


#' @method Nedge multiAphylo_pruner
#' @export
Nedge.multiAphylo_pruner <- Nedge.multiAphylo

#' @method Nnode aphylo
#' @export
Nnode.aphylo <- function(phy, ...) ape::Nnode(phy$tree, ...)

#' @method Nnode aphylo_estimates
#'@export
Nnode.aphylo_estimates <- function(phy, ...) ape::Nnode(phy$dat, ...)

#' @method Nnode multiAphylo
#' @export
Nnode.multiAphylo <- function(phy, ...) sapply(phy, ape::Nnode, ...)

#' @method Nnode aphylo_pruner
#' @export
Nnode.aphylo_pruner <- function(phy, ...) {
  
  dots <- list(...)
  if ("internal.only" %in% names(dots))
    .Nnode_aphylo_pruner(phy, internal_only = dots$internal.only)
  else
    .Nnode_aphylo_pruner(phy, internal_only = TRUE)
  
}

#' @export
#' @method Nnode multiAphylo_pruner
Nnode.multiAphylo_pruner <- Nnode.multiAphylo

#' @export
#' @method Ntip aphylo
Ntip.aphylo <- function(phy) ape::Ntip(phy$tree)

#' @export
#' @method Ntip aphylo_estimates
Ntip.aphylo_estimates <- function(phy) ape::Ntip(phy$dat)

#' @export
#' @method Ntip multiAphylo
Ntip.multiAphylo <- function(phy) sapply(phy, ape::Ntip)

#' @export
#' @method Ntip multiAphylo_pruner
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
#' @examples 
#' # Generating data for the example
#' set.seed(223)
#' dat <- rmultiAphylo(10, n = 5, P = 2)
#' Nann(dat)
#' Nannotated(dat)
#' Ntrees(dat)
NULL

#' @rdname aphylo-info
#' @export
Nann <- function(phy) UseMethod("Nann")

#' @export
# @method Nann aphylo
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
