#' Given that we want to remove a set of \eq{l \subset L}, the algorithm goes as
#' follows:
#' 
#' 1. Match the node ids with position.
#' 2. Sort it to be decreasing, so the most inner nodes show first,
#' 3. Increase the list:
#'   a. Compute the geodesic matrix G,
#'   b. For i in l do, j != i:
#'     If G(i,j)!=0, then add it to the list
#'   
#'   \code{sapply(1:(3-1), function(x) sapply(x:3, c, x)[,-1,drop=FALSE])}
#'   
#'   b. Degine tags(n). For k in p, For s in p[k]:
#'       1. If G[s] != 0, then tag() = 1
#'       
#' Two things to do: 1 remove them from the 
#'       
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
#' plot(prune(x, 5), main="removing 5", show.node.label = TRUE)
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
#' @name prune
NULL

#' @export
#' @rdname prune
prune <- function(x, ids, ...) UseMethod("prune")

#' @export
#' @rdname prune
prune.po_tree <- function(x, ids) {
  
  # 1. Identify which will be removed ------------------------------------------
  
  # Getting the unique set, and sorting it
  ids <- sort(unique(ids))
  n   <- length(attr(x, "labels"))
  
  # Matching to actual labels
  if (is.character(ids))
    ids <- match(ids, getlabels(x)) - 1L
  
  # Checking the lengths
  if (any(is.na(ids)))
    stop("Some -ids- don't match any leafs of the tree.")
  
  if (any(ids > (n - 1)))
    stop("Out of range: Some -ids- are out of range (above n-1).")
  if (any(ids < 1))
    stop("Out of range: Some -ids- are out of range (below 1). Root node cannot be removed.")
  
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
  x <- x[edges_ids,,drop=FALSE]
  
  
  # 4. Relabeling --------------------------------------------------------------
  x[,1] <- new_ids[match(x[,1], old_ids)]
  x[,2] <- new_ids[match(x[,2], old_ids)]
  
  attr(x, "labels") <- old_labels[-(nodes_ids + 1L)]
  
  # 5. Re computing the offspring ----------------------------------------------
  
  attr(x, "offspring") <- list_offspring(x)
  
  structure(x, class= "po_tree")
  
}



