#' Extensions to the `as.phylo` function
#' 
#' This function takes an edgelist and recodes (relabels) the nodes following
#' \CRANpkg{ape}'s coding convention. 
#' 
#' @param x Either an edgelist or an object of class [aphylo].
#' @param ... Further arguments passed to the method.
#' @param edge.length A vector with branch lengths  (optional). 
#' @param root.edge A numeric scalar with the length for the root node (optional).
#' @return An integer matrix of the same dimmension as `edges` with the following
#' aditional attribute:
#' \item{labels}{Named integer vector of size `n`. Original labels of the edgelist
#' where the first `n` are leaf nodes, `n+1` is the root node, and the reminder
#' are the internal nodes.}
#' 
#' @examples 
#' 
#' # A simple example ----------------------------------------------------------
#' # This tree has a coding different from ape's
#' 
#' mytree <- matrix(c(1, 2, 1, 3, 2, 4, 2, 5), byrow = TRUE, ncol=2)
#' mytree
#' 
#' ans <- as.phylo(mytree)
#' ans
#' plot(ans)
#' 
#' 
#' @name as.phylo
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

#' @rdname as.phylo
#' @export
#' @method as.phylo matrix
as.phylo.matrix <- function(
  x,
  edge.length = NULL,
  root.edge   = NULL,
  ...
) {
  
  # Must be a two column matrix
  if (ncol(x) != 2L)
    stop(
      "`x` must be a two-column matrix. It has ", ncol(x), " instead.",
      call. = FALSE
      )
  
  # Dropping incomplete cases
  x <- x[complete.cases(x), , drop = FALSE]
  
  # Retrieving the labels
  if (!inherits(x, "integer")) {
    
    label <- sort(unique(as.vector(x)))
    
    # Recasting as an integer vector
    x <- matrix(match(as.vector(x), label), ncol=2L)
    
  } else {
    
    label <- range(as.vector(x))
    label <- label[1]:label[2]
    
  }
  
  # Computing degrees
  nodes <- sort(unique(as.vector(x)))
  ideg <- tabulate(x[,2] - nodes[1L] + 1L, nbins = nodes[length(nodes)] - nodes[1L] + 1)
  odeg <- tabulate(x[,1] - nodes[1L] + 1L, nbins = nodes[length(nodes)] - nodes[1L] + 1)
  
  # Classifying
  roots <- nodes[ideg == 0 & odeg > 0]
  leafs <- nodes[ideg == 1 & odeg == 0]
  inner <- nodes[ideg == 1 & odeg > 0]
  
  # Multiple parents
  test. <- which(ideg > 1)
  if (length(test.))
    stop("Multiple parents are not supported. The following nodes have multiple parents: ",
         paste(nodes[test.], collapse=", "))
  
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
  test.   <- is.numeric(nodes) &&
    all(sort(leafs) == 1:length(leafs)) &&
    (roots == (nleafs + 1L))
  
  # Case in which the only inner node is the root!
  test. <- ifelse(
    length(inner) > 0,
    test. & (sort(inner) == (nleafs + 2L):length(nodes)),
    test.
    )
    
  # Defining the labels:
  #  - Leafs go from 1 to n
  #  - Root node is n + 1
  #  - And the inner goes from n + 2 to length(nodes)
  # This doest it smoothly
  if (!test.) {
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
    tip.label   = unname(label[leafs]),
    Nnode       = length(inner) + 1L,
    node.label  = unname(label[c(roots, inner)]),
    root.edge   = unname(root.edge)
    )
}


#' @rdname as.phylo
#' @export
#' @method as.phylo aphylo
as.phylo.aphylo <- function(x, ...) {
  
  x$tree
  
}

