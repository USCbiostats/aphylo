#' Recodes an edgelist as a Partially Ordered Tree
#' 
#' The function \code{\link{new_aphylo}} uses this function to make sure that
#' the edgelist that is passed makes a partial order. This is a requirement
#' for the peeling algorithm, which is used explicitly in the \code{LogLike}
#' function.
#' 
#' @details
#' The recoded edgelist is such that in all rows the first element, parent
#' node, as a label that is less than the second element, the offspring, a
#' partial order.
#' 
#' @return A matrix of the same dimension as \code{edges}, an edgelist, recoded
#' to form a partial order. Besides of been of class \code{matrix}, the resulting
#' object is also of class \code{po_tree} and has an aditional attribute:
#' \item{labels}{Named integer vector of size n. Original labels of the edgelist
#' where the names are from 0 to \code{n}.}
#' 
#' @template parameters
#' @templateVar edges 1
#' @export
#' @family Data management functions
#' @examples
#' # Recoding an ape tree -----------------------------------------------------
#' 
#' set.seed(1122233)
#' apetree <- ape::rtree(5)
#' potree  <- as_po_tree(apetree$edge)
#' 
#' apetree$edge
#' potree
#' 
#' # Going back
#' potree[] <- attr(potree, "labels")[potree[] + 1]
#' potree # Ordering is a little off, but is the same tree
#' 
#' @export
as_po_tree <- function(edges) UseMethod("as_po_tree")

#' @export
#' @rdname as_po_tree
as_po_tree.phylo <- function(edges) {
  
  atree   <- edges[["edge"]]

  atree[] <- c(
    edges[["tip.label"]], # 1 ... n
    "root",               # n + 1
    paste0("internal", 1L:(edges[["Nnode"]]))
    )[atree[]]
    
  
  as_po_tree(atree)
}

#' @export
#' @rdname as_po_tree
as_po_tree.default <- function(edges) {
  
  # If PO tree, then nothing to do
  if (inherits(edges, "po_tree"))
    return(edges)
  
  # Checking class
  if (!inherits(edges, c("data.frame", "matrix")))
    stop("-edges- must be either an object of class -data.frame- or -matrix- (is of class ",
      class(edges),")."
    )
  
  # Checking size
  P <- ncol(edges)
  if (P != 2L)
    stop(sprintf("-edges- must have at least 2 columns (it only has %i).", P))
  
  # In the case of data.frame
  if (inherits(edges, "data.frame")) {
    warning("-edges- is of class -data.frame- and will be coerced to a matrix.")
    edges <- as.matrix(edges)
  }
  
  # Checking if it is character
  m <- nrow(edges)
  if (mode(edges) == "character") {
    
    # Retrieving labels
    labels <- as.factor(c(edges[,1,drop=TRUE], edges[,2,drop=TRUE]))
    
    # Getting it as integer
    edges      <- cbind(
      as.integer(labels[1:m]),
      as.integer(labels[(m+1):length(labels)])
    )
    
    # And the corresponding labels
    labels <- levels(labels)
    
  } else {
    
    edges      <- apply(edges, 2, as.integer)
    labels <- NULL
    
  }
  
  # Recoding as a PO tree and recycling the labels
  edges <- recode_as_po(edges)
  if (length(labels)) {
    attr(edges, "labels") <- structure(
      labels[attr(edges, "labels")],
      names = names(attr(edges, "labels"))
    )
  }
  
  edges
  
}

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
#' @family Data management functions
#' @name aphylo-class
NULL

# 
# set.seed(1)
# x <- sim_tree(5)
# x[] <- letters[x[] + 1]
# attr(x, "offspring") <- NULL
# check_edgelist(as.matrix(unclass(x)));x

check_annotations <- function(x) {
  
  # Checking class
  if (!inherits(x, c("data.frame", "matrix")))
    stop("-annotations- must be either an object of class -data.frame- or -matrix- (is of class ",
         class(x),").")
  
  # Checking size
  P <- ncol(x) - 1L
  if (P < 1L)
    stop("-annotations- must have at least 2 columns (it only has ", P + 1L, ").")
  
  # Checking class and coercing
  node_labels <- as.character(x[,1])
  test        <- apply(x[,-1,drop=FALSE], 2, typeof)
  test        <- which(!sapply(test, `%in%`, table = c("integer", "double")))
  if (length(test))
    stop("All annotation columns (but the first one) must be either -integer- or -double-.",
         " The following columns are not: ", paste(test, collapse=", "), ".")
  
  node_annotations <- apply(x[,-1,drop=FALSE], 2, as.integer)
  
  # Checking values
  test <- which(apply(node_annotations, 1, function(y) any(!(y %in% c(0, 1, 9, NA)))))
  if (length(test))
    stop("The following rows of -annotations- have values different from c(0, 1, 9, NA): ",
         paste(test, collapse=", "), ".")
  
  # Replacing NAs
  node_annotations[is.na(node_annotations)] <- 9
  
  # Checking colnames
  fun_names <- colnames(node_annotations)
  
  if (!length(fun_names))
    fun_names <- sprintf("fun%03i", P)
  
  # Returning
  list(
    labels      = unname(node_labels),
    annotations = unname(node_annotations),
    names       = fun_names,
    P           = P
    )
  
}

#' @rdname aphylo-class
#' @export
new_aphylo <- function(
  annotations,
  edges
  ) {
  
  # Basic checks ---------------------------------------------------------------
  
  # Checking annotations
  A <- check_annotations(annotations)
  
  # Checking edgelist
  E <- as_po_tree(edges)
  
  # Are all annotations in the edgelist
  test <- which( !(A[["labels"]] %in% attr(E, "labels")) )
  if (length(test))
    stop("The following -nodes- (annotations) are not in present in -edges-:\n",
         paste(A[["labels"]][test], collapse=", "))
  
  # Filling the gaps in A ------------------------------------------------------
  test  <- which( !(attr(E, "labels") %in% A[["labels"]]))
  added <- NULL
  if (length(test)) {
    
    # Filling the Annotations
    A[["annotations"]] <- rbind(
      A[["annotations"]],
      matrix(NA, nrow = length(test), ncol = ncol(A[["P"]]))
    )
    
    # Labels
    A[["labels"]] <- c(A[["labels"]], attr(E, "labels")[test])
    
    # And declaring which ones where added
    added <- attr(E, "labels")[test]
  }
  
  # Sorting the data
  ord <- match(A[["labels"]], attr(E, "labels"))
  A[["annotations"]] <- A[["annotations"]][ord,]
  
  # Listing offsprings
  O <- list_offspring(E)
  
  # Returning
  as_aphylo(
    annotations = A[["annotations"]],
    fun_names   = A[["names"]], 
    offspring   = O,
    noffspring  = sapply(O, length),
    edges       = E
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
  E <- with(x, as_ape_tree(edges))
  
  # Figuring out labels
  labs <- attr(x, "labels")[match(
    attr(E, "nodes")[which(attr(E, "isleaf"))],
    attr(x[["edges"]], "labels")
    )]
  
  
  structure(list(
    edge        = E,
    edge.length = rep(1, nrow(E)),
    tip.label   = labs,
    Nnode       = sum(!attr(E, "isleaf"))
  ), class = "phylo")
  
}

#' @rdname as.apephylo
#' @export
as.apephylo.po_tree <- function(x, ...) {
  # Cleaning attributes
  attr(x, "noffspring") <- NULL
  attr(x, "offspring")  <- NULL
  
  # Recoding edgelist
  E <- as_ape_tree(x)
  
  # Figuring out labels
  labs <- attr(x, "labels")[match(
    attr(E, "labels")[which(attr(E, "isleaf"))],
    names(attr(x, "labels"))
  )]
  
  structure(list(
    edge        = E,
    edge.length = rep(1, nrow(E)),
    tip.label   = labs,
    Nnode       = sum(!attr(E, "isleaf"))
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
    isleaf = c(rep(TRUE, sum(odeg == 0)), rep(FALSE, length(nodes) - sum(odeg==0))),
    class = c("matrix", "ape_tree")
  )
  
}
