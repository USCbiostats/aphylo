#' List each nodes' offspring or parent
#' 
#' For each node in a tree, the functions `list_offspring` and `list_parents`
#' lists all its offspring and parents, respectively.
#' 
#' @param x An object of class `phylo` or `aphylo`.
#' @return List of length `n` (total number of nodes).
#' 
#' @examples 
#' # A simple example with phylo tree ------------------------------------------
#' 
#' set.seed(4)
#' x <- ape::rtree(10)
#' list_offspring(x)
#' 
#' @export
list_offspring <- function(x) UseMethod("list_offspring")

#' @export
list_offspring.aphylo <- function(x) {
  x$offspring
}

#' @export
list_offspring.phylo <- function(x) {
  o <- split(x$edge[,2], x$edge[,1])
  a <- matrix(list(integer(0L)), ncol=1, nrow = ape::Nnode(x, internal.only = FALSE))
  a[as.integer(names(o))] <- unname(o)
  a[,1]
}

#' @export
#' @rdname list_offspring
list_parents <- function(x) UseMethod("list_parents")

#' @export
list_parents.phylo <- function(x) {
  
  ans <- rbind(x$edge, cbind(list(integer(0)), ape::Ntip(x) + 1))
  ans[order(unlist(ans[,2])),][,1]
  
}

#' @export
list_parents.aphylo <- function(x) {
  
  list_parents(x$tree)
  
}
