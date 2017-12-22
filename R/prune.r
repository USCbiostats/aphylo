#' Removing leafs and nodes from a tree
#' 
#' This function takes one or more nodes/leafs from a given tree and removes them
#' making sure that the position indexes are updated, hence preserving the 
#' [po_tree()] structure.
#' 
#' @param x An object of class `po_tree` or `aphylo`.
#' @param ids Either a vector or a scalar indicating which nodes/leafs to remove.
#' If integer, then its values should be within (0, n-1]. Otherwise, if character
#' it can be the nodes/leafs labels.
#' @param ... Ignored
#' 
#' @return A prunned version of the tree.
#' @family Data management functions
#' @details 
#' 
#' From now we will denote node(s) as either internal node(s) or leaf(s). Given
#' that we want to remove `ids`, the algorithm goes as follows:
#' \enumerate{
#' \item Identifies which nodes are been asked to be removed and
#' checks whether these actually exits.
#' 
#' \item Computes the topological shortest path matrix between pairs of
#' nodes, and if the nodes to be removed are parents of nodes not included
#' in the list, these will be added (the whole branch is removed).
#' 
#' \item Considering which nodes are been removed, a new set of positions
#' is computed so that it follows [po_tree()] convention.
#' 
#' \item The edgelist is updated, as well the labels.
#' 
#' \item The set of offspring is recalculated.
#' }
#' 
#' @examples 
#' 
#' # A simple example of how prune works-------------------------------------------
#' # Simulating a nice tree
#' set.seed(1213)
#' x <- sim_tree(4)
#' 
#' # Setting up the plot envir
#' oldpar <- par(no.readonly = TRUE)
#' par(mfrow = c(3,2), mai=rep(.5,4), cex=1, xpd=NA, omi=c(0,0,1,0))
#' 
#' # Plotting 
#' plot(x, main = "Full tree", show.node.label = TRUE)
#' plot(prune(x, c(2,6)), main="removing (2,6)", show.node.label = TRUE)
#' plot(prune(x, 6), main="removing 6", show.node.label = TRUE)
#' plot(prune(x, 4), main="removing 4", show.node.label = TRUE)
#' plot(prune(x, 3), main="removing 3", show.node.label = TRUE)
#' plot(prune(x, c(4,6,3)), main="removing (4,6,3)", show.node.label = TRUE)
#' 
#' # Adding a title
#' par(mai=rep(1,4), omi=rep(0,4), mfrow=c(1,1), new=FALSE)
#' title("Prunning trees with -prune-")
#' par(oldpar)
#' 
#' # Removing the leafs --------------------------------------------------------
#' 
#' set.seed(1)
#' x <- sim_tree(25)
#' oldpar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2))
#' plot(x)
#' plot(prune(x, "leafs"))
#' plot(prune(prune(x, "leafs"), "leafs"))
#' plot(prune(prune(prune(x, "leafs"), "leafs"), "leafs"))
#' par(oldpar)
#' 
#' @name prune
NULL

#' @export
#' @rdname prune
prune <- function(x, ids, ...) UseMethod("prune")

#' @export
#' @rdname prune
prune.po_tree <- function(x, ids, ...) {
  
  # 1. Identify which will be removed ------------------------------------------
  
  # Getting the unique set, and sorting it
  n   <- length(attr(x, "labels"))
  
  # Matching to actual labels
  map_ids_to_positions.po_tree("ids", "x")
  
  if (length(ids) == 0) {
    warning("Nothing to prune.")
    return(x)
  }
  
  # 2. Computing Geodesics to extend the list ----------------------------------
  nodes_ids <- ids
  G <- approx_geodesic(x, undirected = FALSE, warn = FALSE)
  
  # Which others should be removed
  for (l in ids)
    nodes_ids <- c(nodes_ids, which(G[l + 1L, ] > 0) - 1L)
  
  # Getting the final list and checking if it leaves at least 2 nodes
  nodes_ids <- sort(unique(nodes_ids))
  
  if ( (n - length(nodes_ids)) < 2 )
    stop("You are removing all but the root node, and that's not a tree.")
  
  # 3. Marking the ones that need to be removed --------------------------------
  old_ids <- 0L:(n - 1L)
  new_ids <- rep(0L, n)
  new_ids[nodes_ids+1L] <- 1L
  new_ids <- old_ids - cumsum(new_ids)
  
  old_labels <- attr(x, "labels")
  
  # 4. Removing everything that won't be used any more -------------------------
  
  # From the edgelist
  edges_ids <- which(!(x[,1] %in% nodes_ids) & !(x[,2] %in% nodes_ids))
  edge.length <- attr(x, "edge.length")[edges_ids]
  x <- x[edges_ids,,drop=FALSE]

  # 4. Relabeling --------------------------------------------------------------
  x[,1] <- new_ids[match(x[,1], old_ids)]
  x[,2] <- new_ids[match(x[,2], old_ids)]
  
  names(old_labels) <- new_ids

  # 5. Re computing the offspring ----------------------------------------------
  new_po_tree(
    edges       = x,
    edge.length = edge.length,
    labels      = old_labels[-(nodes_ids + 1L)]
    )
  
}

#' @rdname prune
#' @export
prune.phylo <- function(x, ids, ...) {
  as.apephylo(prune(as_po_tree(x), ids))
}

