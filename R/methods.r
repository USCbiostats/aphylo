#' Plot and print methods for \code{phylo_offspring} objects
#' @param x An object of class \code{phylo_offspring}.
#' @param y Ignored.
#' @param tip.color See \code{\link[ape:plot.phylo]{plot.phylo}}
#' @param ... Further arguments passed to the method.
#' @rdname get_offspring
#' @export
plot.phylo_offspring <- function(
  x, y=NULL, tip.color=NULL, ...) {
  
  if (!length(tip.color)) {
    tip.color <- colors(9)[with(x, experiment[leaf_node,1,drop=TRUE])]
  } 
  
  plot(as.phylo.phylo_offspring(x), tip.color=tip.color, ...)
  
}