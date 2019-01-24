#' APE wrappers
#' @param phy An object of class [aphylo] (see [ape::Nnode]).
#' @param ... Further arguments passed to the method.
#' 
#' Functions `Nedge`, `Nnode`, and `Ntip` have S3 methods in the `ape` package.
#' 
#' The function `Nann`, which reports the number of annotations that the tree has,
#' is not a function present in ape, but has a method for
#' both `phylo` and `multiPhylo` objects which returns a warning an 0L.
#' @return Integer with the number of edges, nodes, or tips accordignly.
#' @name ape-methods
NULL

#' @rdname ape-methods
#' @export
Nedge.aphylo <- function(phy) ape::Nedge(phy$tree)

#' @rdname ape-methods
#' @export
Nnode.aphylo <- function(phy, ...) ape::Nnode(phy$tree, ...)

#' @rdname ape-methods
#' @export
Ntip.aphylo <- function(phy) ape::Ntip(phy$tree)



#' @rdname ape-methods
#' @export
Nann <- function(phy) UseMethod("Nann")

#' @rdname ape-methods
#' @export
Nann.aphylo <- function(phy) ncol(phy$tip.annotation)

#' @rdname ape-methods
#' @export
Nann.phylo <- function(phy) {
  warning("phylo objects have no annotations.", call. = FALSE)
  0L
}

#' @rdname ape-methods
#' @export
Nann.multiPhylo <- function(phy) {
  warning("multiPhylo objects have no annotations.", call. = FALSE)
  0L
}



#' Indexing aphylo objects
#' @param x An object of class [aphylo].
#' @param i,j Integer vector. Indices of genes or functions.
#' @param value Integer vector. Replacing values, can be either `c(0, 1, 9, NA)`.
#' @details The subsetting method allows selecting one or more annotations from
#' the [aphylo] object.
#' @name aphylo-index
NULL

#' @export
#' @rdname aphylo-index
`[.aphylo` <- function(x, j) {
  

  x$tip.annotation  <- x$tip.annotation[, j, drop=FALSE]
  x$node.annotation <- x$node.annotation[, j, drop=FALSE]
  
  x
  
}


#' @export
#' @rdname aphylo-index
`[<-.aphylo` <- function(x, i, j, value) {
  
  # Checking all values
  value[is.na(value)] <- 9L
  test <- which(!(value %in% c(0L, 1L, 9L)))
  if (length(test)) 
    stop(
      "Replacing annotations must be either c(0, 1, 9, NA).",
      " The following values are not valid: ",
      paste(value[test], collapse=", "), ".", call. = FALSE
      )
  
  if (missing(j))
    j <- 1L:Nann(x)
  if (missing(i))
    i <- 1L:Nnode(x, internal.only = FALSE)
  
  # Sorting
  nt <- Ntip(x)
  tips_ids  <- which(i <= nt)
  nodes_ids <- which(i > nt)
  
  if (length(tips_ids))
    x$tip.annotation[i[tips_ids], j] <- value[tips_ids]
  
  if (length(nodes_ids))
    x$node.annotation[i[nodes_ids] - nt, j] <- value[nodes_ids]
  
  x
}