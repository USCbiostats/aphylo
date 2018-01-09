#' Extensions to the `as.phylo` function
#' 
#' This function takes an edgelist and recodes (relabels) the nodes following
#' \CRANpkg{ape}'s coding convention. `as_ape_tree` is the powerhorse of
#' [as.apephylo()].
#' 
#' @template parameters
#' @templateVar edges 1
#' @return An integer matrix of the same dimmension as `edges` with the following
#' aditional attribute:
#' \item{labels}{Named integer vector of size `n`. Original labels of the edgelist
#' where the first `n` are leaf nodes, `n+1` is the root node, and the reminder
#' are the internal nodes.}
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
#' @name as.phylo.matrix
NULL

#' Creates a phylo object
#' @noRd
new_phylo <- function(
  edge,
  tip.label,
  Nnode,
  edge.length = NULL,
  node.label  = NULL,
  root.edge   = NULL
) {
  
  structure(
    c(
      list(edge = edge),
      # Since edge.length is optional
      if (length(edge.length))
        list(edge.length = edge.length)      
      else 
        NULL,
      list(
        tip.label  = tip.label,
        Nnode      = Nnode,
        node.label = node.label
      ),
      # Since root.edge is optional
      if (length(root.edge))
        list(root.edge = root.edge)
      else
        NULL      
    ),
    class = "phylo"
  )
}


#' @rdname as.phylo.matrix
#' @export
as.phylo.matrix <- function(
  x,
  edge.length = NULL,
  root.edge   = NULL,
  ...
  ) {
  
  # Computing degrees
  nodes <- unique(as.vector(x))
  ideg  <- fast_table_using_labels(x[,2], nodes)
  odeg  <- fast_table_using_labels(x[,1], nodes)
  
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
  
  # We will not relabel iff:
  # 1. nodes is integer/numeric vector
  # 2. Leafs are continuously labeled from 1 to n
  # 3. Root is n+1
  # 4. Interior nodes are from n+2 to m
  nleafs <- length(leafs)
  test   <- is.numeric(nodes) &&
    all(sort(leafs) == 1:length(leafs)) &&
    (roots == (nleafs + 1L)) &&
    (sort(inner) == (nleafs + 2L):length(nodes))
  
  # Defining the labels:
  #  - Leafs go from 1 to n
  #  - Root node is n + 1
  #  - And the inner goes from n + 2 to length(nodes)
  # This doest it smoothly
  if (!test) {
    nodes <- c(leafs, roots, inner)
    
    # Finding indexes and corresponding new labels
    iroots <- which(x[] == roots)
    lroots <- match(roots, nodes)
    
    ileafs <- which(x[] %in% leafs)
    lleafs <- match(x[ileafs], nodes)
    
    iinner <- which(x[] %in% inner)
    linner <- match(x[iinner], nodes)
    
    # Relabeling the edgelist
    x[iroots] <- lroots
    x[ileafs] <- lleafs
    x[iinner] <- linner
  }
  
  # Returning the `phylo` object
  new_phylo(
    edge        = unname(x),
    edge.length = unname(edge.length),
    tip.label   = unname(leafs),
    Nnode       = length(inner) + 1L,
    node.label  = unname(c(roots, inner)),
    root.edge   = unname(root.edge)
    )
}


#' Coercing into `phylo` objects of the \pkg{ape} package.
#'
#' `as.apephylo` coerces objects to [ape:as.phylo::phylo()] objects
#' from the \pkg{ape} package. These have several methods including visualization
#' methods that can be useful.
#'
#' @param x An object of class [=as_po_tree::po_tree()] or
#' [=new_aphylo::aphylo()].
#' @param ... Ignored.
#' @return An object of class [ape:as.phylo::phylo()]
#' @family Data management functions
#' @export
as.phylo.aphylo <- function(x, ...) {
  
  x$tree
  
}

