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
#' as_po_tree(apetree)
#' 
#' @export
as_po_tree <- function(edges) UseMethod("as_po_tree")

#' @export
#' @rdname as_po_tree
as_po_tree.phylo <- function(edges) {
  
  atree    <- edges[["edge"]]

  nodelabs <- edges[["node.label"]]
  if (!length(nodelabs))
    nodelabs <- c("root", paste0("internal", 1L:(edges[["Nnode"]])))
  
  atree[] <- c(
    edges[["tip.label"]], # 1 ... n
    nodelabs
    )[atree[]]
    
  
  as_po_tree(atree)
}

#' @rdname as_po_tree
#' @export
#' @return In the case of \code{is.po_tree}, \code{TRUE} if is an object of
#' class \code{po_tree}, \code{FALSE} otherwise.
is.po_tree <- function(x) inherits(x, "po_tree")

#' @export
#' @rdname as_po_tree
print.po_tree <- function(x, ...) {
  # Which are leafs
  dgr <- fast_table_using_labels(
    x[,1],
    as.integer(names(attr(x, "labels")))
    )
  
  cat(
    "\nA PARTIALLY ORDERED PHYLOGENETIC TREE\n",
    sprintf("  # Internal nodes: %i", attr(x, "Nnode")),
    sprintf("  # Leaf nodes    : %i\n", length(attr(x, "labels")) - attr(x, "Nnode")),
    sprintf("  Leaf nodes labels: \n    %s, ...\n",
            paste0(utils::head(attr(x, "labels")[which(dgr==0)]), collapse = ", ")
            ),
    sprintf("  Internal nodes labels:\n    %s, ...\n",
            paste0(utils::head(attr(x, "labels")[which(dgr>0)]), collapse = ", ")
    ),
    
    sep = "\n"
  )
  
  invisible(x)
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
  
  # Treating as character
  m <- nrow(edges)
  edges <- apply(edges, 2, as.character)

  # Retrieving labels
  labels <- as.factor(c(edges[,1,drop=TRUE], edges[,2,drop=TRUE]))
  
  # Getting it as integer
  edges      <- matrix(as.integer(labels), byrow = FALSE, ncol = 2)
  
  # And the corresponding labels
  labels <- levels(labels)
  
  # Recoding as a PO tree and recycling the labels
  edges <- recode_as_po(edges)
  if (length(labels)) {
    attr(edges, "labels") <- structure(
      labels[as.integer(attr(edges, "labels"))],
      names = names(attr(edges, "labels"))
    )
  }
  
  edges
  
}

#' Annotated Phylogenetic Tree
#' 
#' The \code{aphylo} class tree holds both the tree structure represented as a
#' partially ordered phylogenetic tree, and node annotations. While annotations
#' are included for both leafs and inner nodes, the algorithms included in this
#' package only uses the leaf annotations.
#' 
#' @template parameters
#' @templateVar edges 1
#' @templateVar annotationslab 1
#' @param ... Further parameters passed to the method.
#' 
#' @details Plotting is done via \code{\link[ape:plot.phylo]{plot.phylo}} 
#' from the \CRANpkg{ape} package.
#' 
#' @return A list of class \code{aphylo} with the following elements:
#' \item{annotations}{A matrix of size \eqn{n\times P}{n * P} with leaf functional
#' annotations. These can be either 0, 1, or 9 indicating no function, function, and
#' no information respectively.}
#' \item{offspring}{A list of size \eqn{n}. The \eqn{i}-th element of the list
#' can be either \code{NULL} (meaning no offspring), or a vector of length \code{noffspring[i]},
#' with the relative positions of the offspring from \eqn{1} to \eqn{n}.} 
#' \item{noffspring}{An integer vector of size \eqn{n}. Indicates the number of
#' offspring each node has.}
#' \item{edges}{An integer matrix of class \code{po_tree}. An edgelist (see 
#' \code{\link{as_po_tree}}).}
#' 
#' @examples 
#' # A simple example ----------------------------------------------------------
#' 
#' data(fakeexperiment)
#' data(faketree)
#' ans <- new_aphylo(fakeexperiment, faketree)
#'  
#' # We can visualize it
#' plot(ans)
#' @family Data management functions
#' @family aphylo methods
#' @name aphylo-class
NULL

# 
# set.seed(1)
# x <- sim_tree(5)
# x[] <- letters[x[] + 1]
# attr(x, "offspring") <- NULL
# check_edgelist(as.matrix(unclass(x)));x

duplicate_check <- function(x, what) {
  f <- as.factor(x)
  l <- levels(f)
  d <- fast_table(as.integer(f))
  ans <- l[d[which(d[,2] > 1),1]]
  if (length(ans))
    stop("The following ", what, " have repeated measures:\n  ",
         paste(sprintf("- %s (x%i)", ans, d[which(d[,2] > 1),2]), collapse=",\n  "), ".")
  
  invisible(0)
}

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
  duplicate_check(node_labels, "genes in -annotations-")
  
  
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
  
  # Setting dimnames
  dimnames(node_annotations) <- list(
    unname(node_labels),
    fun_names
  )
  
  # Actual annotations
  node_annotations
  
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
  test <- which( !(rownames(A) %in% attr(E, "labels")) )
  if (length(test))
    stop("The following -nodes- (annotations) are not in present in -edges-:\n",
         paste(rownames(A)[test], collapse=", "))
  
  # Filling the gaps in A ------------------------------------------------------
  test  <- which( !(attr(E, "labels") %in% rownames(A)))
  added <- NULL
  if (length(test)) {
    
    # Filling the Annotations
    A <- rbind(
      A,
      matrix(9, nrow = length(test), ncol = ncol(A),
             dimnames = list(
               attr(E, "labels")[test],
               colnames(A)
             ))
    )
    
    # And declaring which ones where added
    added <- attr(E, "labels")[test]
  }
  
  # Sorting the data
  A <- A[attr(E, "labels"), ,drop=FALSE]
  
  # Listing offsprings
  O <- list_offspring(E)
  
  # Returning
  as_aphylo(
    annotations = A,
    offspring   = O,
    noffspring  = sapply(O, length),
    edges       = E
  )
  
}

as_aphylo <- function(
  annotations,
  offspring,
  noffspring,
  edges,
  checks     = TRUE
  ) {
  
  if (checks) {
    
    # Checking class
    stopifnot(is.matrix(annotations))
    stopifnot(is.list(offspring))
    stopifnot(is.vector(noffspring))
    stopifnot(is.po_tree(edges))
    
    # Checking dimmensions
    d <- dim(annotations)
    stopifnot(d[1] == length(offspring))
    stopifnot(d[1] == length(noffspring))
    stopifnot(length(unique(as.vector(edges))) == d[1])
    
  }
  
  structure(
    list(
      annotations = annotations,
      offspring   = offspring,
      noffspring  = noffspring,
      edges       = edges
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
  tiplabs <- attr(x[["edges"]], "labels")[match(
    attr(E, "labels")[which(attr(E, "isleaf"))],
    attr(x[["edges"]], "labels")
    )]
  
  # Internal nodes labels
  nodelabs <- attr(x[["edges"]], "labels")[match(
    attr(E, "labels")[which(!attr(E, "isleaf"))],
    attr(x[["edges"]], "labels")
  )]
  
  
  structure(list(
    edge        = E,
    edge.length = rep(1, nrow(E)),
    tip.label   = as.character(unname(tiplabs)),
    Nnode       = sum(!attr(E, "isleaf")),
    node.label  = as.character(unname(nodelabs))
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
  tiplabs <- attr(x, "labels")[match(
    attr(E, "labels")[which(attr(E, "isleaf"))],
    names(attr(x, "labels"))
  )]
  
  nodelabs <- attr(x, "labels")[match(
    attr(E, "labels")[which(!attr(E, "isleaf"))],
    names(attr(x, "labels"))
  )]
  
  structure(list(
    edge        = E,
    edge.length = rep(1, nrow(E)),
    tip.label   = as.character(unname(tiplabs)),
    node.label  = as.character(unname(nodelabs)),
    Nnode       = sum(!attr(E, "isleaf"))
  ), class = "phylo")
}


#' Plot and print methods for \code{aphylo} objects
#' 
#' @param x An object of class \code{aphylo}.
#' @param y Ignored.
#' @param geom.tiplab.args Further arguments passed to \code{\link[ggtree:ggtree]{ggtree}}
#' @param gheatmap.args Further arguments passed to \code{\link[ggtree:ggtree]{ggtree}}
#' @param scale.fill.args Further arguments passed to \code{\link[ggtree:ggtree]{ggtree}}
#' @param ... Further arguments passed to the method.
#' @name aphylo-methods
#' @details The \code{plot.aphylo} function is a wrapper of \code{\link[ggtree:ggtree]{ggtree}}
#' that creates a visualization as follows:
#' \enumerate{
#'   \item Retrieve the annotations from the \code{\link[=aphylo-class]{aphylo}} object
#'   \item Create a \code{ggtree} map adding \code{\link[ggtree:geom_tiplab]{geom_tiplab}}s
#'   \item Use \code{\link[ggtree:gheatmap]{gheatmap}} to add a heatmap.
#'   \item Set the colors using \code{\link[ggplot2:scale_fill_manual]{scale_fill_manual}}
#' }
#' 
#' @export
#' @return In the case of \code{plot.aphylo}, an object of class \code{c("ggtree", "gg", "ggplot")}
#' @family aphylo methods
plot.aphylo <- function(
  x,
  y=NULL,
  geom.tiplab.args = list(),
  gheatmap.args    = list(),
  scale.fill.args  = list(), 
  ...
) {
  
  # Checkingout arguments
  if (!("align" %in% names(geom.tiplab.args))) geom.tiplab.args$align <- TRUE
  
  if (!("width" %in% names(gheatmap.args))) gheatmap.args$width <- .25
  if (!("color" %in% names(gheatmap.args))) gheatmap.args$color <- "transparent"
  
  if (!("breaks" %in% names(scale.fill.args))) scale.fill.args$breaks <- c("0", "1", "9")
  if (!("values" %in% names(scale.fill.args))) scale.fill.args$values <- c("gray", "steelblue", "white")
  if (!("labels" %in% names(scale.fill.args))) scale.fill.args$labels <- c("No function", "Function", "N/A")
  
  # Retrieving the annotations
  A <- as.data.frame(apply(x$annotations,2,as.character))
  
  # Creating the mapping
  p <- ggtree::ggtree(as.apephylo(x$edges)) +
    do.call(ggtree::geom_tiplab, geom.tiplab.args)
  
  # Adding the functions
  p <- do.call(ggtree::gheatmap, c(list(p=p, data=A),gheatmap.args)) +
    do.call(ggplot2::scale_fill_manual, scale.fill.args)
  
  p
}

#' @rdname aphylo-methods
#' @export
leafs <- function(x) {
  if (!inherits(x, "aphylo"))
    stop("No -leafs- method for class -", class(x),"-.")
  
  labs <- attr(x$edges, "labels")
  d <- fast_table_using_labels(x$edges[,1], as.integer(names(labs)))
  
  unname(labs[which(d == 0)])
}

#' @export
#' @rdname aphylo-methods
#' @return In the case of \code{print.aphylo}, the same \code{aphylo} object
#' invisible.
print.aphylo <- function(x, ...) {
  
  print(x$edges)
  
  cat("ANNOTATIONS:\n")
  
  ids <- leafs(x)
  print(utils::head(x$annotations[ids,,drop=FALSE]))
  
  if (length(ids) > 6)
    cat("\n...(", length(ids) - 6, " obs. omitted)...\n\n", sep="")
  
  
  invisible(x)
}

#' @export
#' @param object An object of class \code{aphylo}.
#' @rdname aphylo-methods
summary.aphylo <- function(object, ...) {
  ans <- apply(object$annotations, 2, table)
  
  if (!inherits(ans, "list"))
    ans <- list(ans[,1,drop=TRUE])
  
  ans <- do.call(
    rbind, 
    lapply(ans, function(x)
    data.frame(
      `0` = x["0"],
      `1` = x["1"],
      `NA` = x["9"],
      check.names = FALSE
      )
    )
  )
  rownames(ans) <- colnames(object$annotations)
  cat("\nDistribution of functions:\n")
  print(ans)
  
  invisible(ans)
  
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
#' edges <- sim_tree(30)$edges
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
