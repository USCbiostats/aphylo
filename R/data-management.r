#' Annotated Phylogenetic Tree
#' 
#' 
#' 
#' @param leafidvar A character scalar with the name of the leaf id variable
#' in \code{annotations}.
#' @template parameters
#' @templateVar edges 1
#' @templateVar annotations 1
#' @param funvars A character vector with the names of the function indicator
#' (0-1) variables in \code{data_exper}. If not provided, it is assumed that
#' all but \code{leafidvar} are function variables.
#' 
#' @details Plotting is done via \code{\link[ape:plot.phylo]{plot.phylo}} 
#' from the \CRANpkg{ape} package.
#' 
#' @return A list of length \eqn{n} with relative position of offspring
#' of each node, with respect to \code{data_exper}, starting from 0.
#' @examples 
#' # A simple example ----------------------------------------------------------
#' 
#' data(fakeexperiment)
#' data(faketree)
#' ans <- new_aphylo(fakeexperiment, faketree, "LeafId")
#'  
#' # We can visualize it
#' plot(ans)
#' @export
#' @family Data management functions
new_aphylo <- function(
  annotations,
  edges,
  leafidvar,
  funvars=NULL
  ) {
  
  # Step 0: Do the variables exists?
  if (!(leafidvar %in% colnames(annotations)))
    stop("The variable -", leafidvar, "- does not exists in -annotations-.")
  
  # First check: All nodes in annotations must be in the tree!
  test <- which(!(annotations[[leafidvar]] %in% as.vector(edges)))
  if (length(test)) {
    stop("The following nodes (",length(test),"/",nrow(annotations),
         ") do not show in -edges- (see ?new_aphylo):\n",
         paste(annotations[[leafidvar]][test], collapse = ", "))
  }
  
  # Identifying function variables
  if (!length(funvars))
    funvars <- colnames(annotations)[colnames(annotations) != leafidvar]
  else {
    test <- which(!(funvars %in% colnames(annotations)))
    if (length(test))
      stop("The following elements of -funvars- are not in -annotations-:\n",
           paste(funvars[test], collapse = ", "))
  }
  
  # Coercing edges to aphylo
  if (!inherits(edges, "po_tree")) 
    aphylo <- as_po_tree(edges)
  else 
    aphylo <- edges
  
  labels <- attr(aphylo, "labels")
  
  # Extending and sorting annotations
  dat <- data.frame(pos = 0L:(length(labels) - 1L), id = labels)
  dat <- merge(dat, annotations, by.x="id", by.y=leafidvar, all=TRUE)
  dat <- subset(
    dat[order(dat$pos),],
    select=colnames(annotations)[colnames(annotations) != leafidvar])
  
  dat <- as.matrix(dat)
  dat[is.na(dat)] <- 9
  
  # Listing offsprings
  ans <- list_offspring(aphylo)
  
  # Returning
  as_aphylo(
    annotations = dat,
    fun_names  = funvars, 
    offspring  = ans,
    noffspring = sapply(ans, length),
    edges      = aphylo
  )
  
}

as_aphylo <- function(
  annotations,
  fun_names,
  offspring,
  noffspring,
  edges,
  checks = FALSE) {
  
  structure(
    list(
      annotations = annotations,
      fun_names  = fun_names,
      offspring  = offspring,
      noffspring = noffspring,
      edges      = edges
    ),
    class = "aphylo"
  )
}

#' Coercing into \code{phylo} objects of the \pkg{ape} package.
#'
#' \code{as.apephylo} coerces objects to \code{\link[ape:as.phylo]{phylo}} objects
#' from the \pkg{ape} package. These have several methods including visualization
#' methods that can be useful.
#'
#' @param x An object of class \code{\link[=as_po_tree]{po_tree}} or
#' \code{\link[=new_aphylo]{aphylo}}.
#' @param ... Ignored.
#' @return An object of class \code{\link[ape:as.phylo]{phylo}}
#' @family Data management functions
#' @export
#' @aliases as.aphylo
as.apephylo <- function(x, ...) UseMethod("as.apephylo")

#' @rdname as.apephylo
#' @export
as.apephylo.aphylo <- function(x, ...) {
  
  # Recoding edgelist
  x <- with(x, as_ape_tree(edges))
  
  structure(list(
    edge        = x,
    edge.length = rep(1, nrow(x)),
    tip.label   = sprintf("leaf%03i", 1:sum(attr(x, "isleaf"))),
    Nnode       = sum(!attr(x, "isleaf"))
  ), class = "phylo")
  
}

#' @rdname as.apephylo
#' @export
as.apephylo.po_tree <- function(x, ...) {
  # Cleaning attributes
  attr(x, "noffspring") <- NULL
  attr(x, "offspring")  <- NULL
  
  x <- as_ape_tree(x)
  
  structure(list(
    edge        = x,
    edge.length = rep(1, nrow(x)),
    tip.label   = sprintf("leaf%03i", 1:sum(attr(x, "isleaf"))),
    Nnode       = sum(!attr(x, "isleaf"))
  ), class = "phylo")
}


#' Plot and print methods for \code{aphylo} objects
#' @param x An object of class \code{aphylo}.
#' @param y Ignored.
#' @param tip.color See \code{\link[ape:plot.phylo]{plot.phylo}}
#' @param ... Further arguments passed to the method.
#' @rdname new_aphylo
#' @export
plot.aphylo <- function(
  x, y=NULL, tip.color=NULL, ...) {
  # 
  # if (!length(tip.color)) {
  #   tip.color <- colors(9)[with(x, experiment[leaf_node,1,drop=TRUE])]
  # } 
  # 
  plot(as.apephylo.aphylo(x), tip.color=tip.color, ...)
  
}

#' @param x An object of class \code{aphylo}.
#' @param y Ignored.
#' @param tip.color See \code{\link[ape:plot.phylo]{plot.phylo}}
#' @param ... Further arguments passed to the method.
#' @rdname as_po_tree
#' @export
plot.po_tree <- function(
  x, y=NULL, tip.color=NULL, ...) {

  plot(as.apephylo(x), tip.color=tip.color, ...)
  
}

#' Relabel a tree (edgelist) using \pkg{ape}'s convention
#' 
#' \code{as_ape_tree} is the powerhorse of \code{\link{as.apephylo}}.
#' 
#' @template parameters
#' @templateVar edges 1
#' @return An integer matrix of the same dimmension as \code{edges} with the following
#' two aditional attributes:
#' \item{labels}{Named integer vector of size \code{n}. Original labels of the edgelist
#' where the first \code{n} are leaf nodes, \code{n+1} is the root node, and the reminder
#' are the internal nodes.}
#' \item{isleaf}{Logical vector of size \code{n}. \code{TRUE} when the node is a leaf node.}
#' 
#' @examples 
#' 
#' # A simple example ----------------------------------------------------------
#' set.seed(1231)
#' 
#' # Generating a random tree: Observe that the nodes are coded such that
#' # for all (i,j) that are offspring and parent i > j. This is a requirement
#' # of the peeling algorithm
#' edges <- sim_tree(30)
#' 
#' # In the ape package, nodes are labeled such that leafs are from 1 to n
#' # root node is n+1 and interior nodes can have whatever label.
#' ape_edges <- as_ape_tree(edges)
#' 
#' # And further, is easy to go back to the original code
#' edges_org <- ape_edges
#' edges_org[] <- attr(ape_edges, "labels")[ape_edges[]]
#' 
#' # Comparing the three versions
#' head(edges)
#' head(ape_edges)
#' head(edges_org)
#' 
#' @export
#' @rdname as.apephylo
as_ape_tree <- function(edges) {
  
  # Computing degrees
  nodes <- unique(as.vector(edges))
  ideg  <- fast_table_using_labels(edges[,2], nodes)
  odeg  <- fast_table_using_labels(edges[,1], nodes)
  
  # Classifying
  roots <- nodes[ideg == 0 & odeg > 0]
  leafs <- nodes[ideg == 1 & odeg == 0]
  inner <- nodes[ideg == 1 & odeg > 0]
  
  # Multiple parents
  test <- which(ideg > 1)
  if (length(test))
    stop("Multiple parents are not supported. The following nodes have multiple parents: ",
         paste(nodes[test], collapse=", "))

  # Finding roots
  if (length(roots) > 1)
    stop("Multiple root nodes are not supported.")
  if (length(roots) == 0)
    stop("Can't find a root node here.")
  
  # Defining the labels:
  #  - Leafs go from 1 to n
  #  - Root node is n + 1
  #  - And the inner goes from n + 2 to length(nodes)
  # This doest it smoothly
  nodes <- structure(c(leafs, roots, inner), names = 1:length(nodes))
  
  # Finding indexes and corresponding new labels
  iroots <- which(edges[] == roots)
  lroots <- match(roots, nodes)
  
  ileafs <- which(edges[] %in% leafs)
  lleafs <- match(edges[ileafs], nodes)
  
  iinner <- which(edges[] %in% inner)
  linner <- match(edges[iinner], nodes)
  
  # Relabeling the edgelist
  edges[iroots] <- lroots
  edges[ileafs] <- lleafs
  edges[iinner] <- linner
  
  structure(
    edges,
    labels = nodes,
    isleaf = odeg == 0,
    class = c("matrix", "ape_tree")
  )
  
}
