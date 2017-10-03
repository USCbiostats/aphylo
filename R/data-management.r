# This function takes and edgelist and tries to find 'labels' attr on it. If
# there are no labels, it will make some.
# This only work for matrices (internal use only)
getlabels <- function(edgelist) {
  
  # Checking the class
  if (!is.integer(edgelist))
    stop("-edgelist- must be of class integer.")
  
  # Are there any labels?
  labels <- attr(edgelist, "labels")
  if (length(labels))
    return(unname(labels))
  
  warning("No labels found. We will use the coding as labels.")
  
  # If there are no labels, then we have to 'find them'
  sort(unique(as.vector(edgelist)))
}


# Returns a logical vector indicating whether the node is a leaf or not
# The edgelist must be codded from 0 to n-1.
isleaf <- function(edgelist, from0=TRUE) {
  
  n <- length(getlabels(edgelist))
  
  # Tabulating the parents, if returns 0 means that it has no
  # offspring, hence, is a leaf.
  fast_table_using_labels(edgelist[,1], (1L-from0):(n - from0)) == 0
}


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
#' object is also of class \code{po_tree} and has an aditional attributes:
#' \item{Nnode}{Integer scalar. Number of leaf nodes.}
#' \item{edge.length}{Numeric vector of length \code{nrow(edges)}. Branch lengths.}
#' \item{labels}{Character vector of size n. Original labels of the edgelist.}
#' \item{offspring}{A list of length \eqn{n} with integer vectors. List the
#' offsprings that each node has relative to 0, where 0 is the root node.}
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
#' potree  <- as_po_tree(apetree)
#' 
#' potree
#' 
#' # Going back
#' as.apephylo(potree)
#'
#' @aliases po_tree 
#' @export
as_po_tree <- function(edges, ...) UseMethod("as_po_tree")

#' @export
#' @rdname as_po_tree
#' @details \code{as.po_tree} is an alias for \code{as_po_tree}
as.po_tree <- as_po_tree

#' @export
#' @rdname as_po_tree
as_po_tree.phylo <- function(edges, ...) {
  
  atree    <- edges[["edge"]]

  nodelabs <- edges[["node.label"]]
  if (!length(nodelabs))
    nodelabs <- c("root", paste0("internal", 1L:(edges[["Nnode"]])))
  
  atree[] <- c(
    edges[["tip.label"]], # 1 ... n
    nodelabs
    )[atree[]]
    
  
  as_po_tree(atree, edge.length = edges[["edge.length"]])
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
    x[,1L],
    0L:(length(attr(x, "labels")) - 1L)
    )
  
  leafs <- dgr == 0
  
  cat(
    "\nA PARTIALLY ORDERED PHYLOGENETIC TREE\n",
    sprintf("  # Internal nodes: %i", attr(x, "Nnode")),
    sprintf("  # Leaf nodes    : %i\n", length(attr(x, "labels")) - attr(x, "Nnode")),
    sprintf("  Leaf nodes labels: \n    %s%s\n",
            paste0(utils::head(attr(x, "labels")[which(leafs)]), collapse = ", "),
            ifelse(sum(leafs) > 6, ", ...", ".")
            ),
    sprintf("  Internal nodes labels:\n    %s%s\n",
            paste0(utils::head(attr(x, "labels")[which(!leafs)]), collapse = ", "),
            ifelse(sum(!leafs) > 6, ", ...", ".")
    ),
    
    sep = "\n"
  )
  
  invisible(x)
}

# This function is intended to be used internally
new_po_tree <- function(
  edges,
  labels, 
  edge.length        = NULL,
  offspring       = NULL,
  check.edges     = TRUE,
  check.edge.length  = TRUE,
  check.labels    = TRUE,
  check.offspring = TRUE
  ) {
  
  # Edges ----------------------------------------------------------------------
  
  # Since when checking edges we are able to id -n-, we have to do it all
  # the time that some other thing that needs is checked.
  
  if (check.edges | check.edge.length | check.labels | check.offspring) {
    
    # Should be a matrix
    if (!is.matrix(edges))
      stop("-edges- should be an object of class matrix")
    
    # Coercion into integers
    if (!is.integer(edges))
      edges[] <- as.integer(edges[])
    
    # Are there any NAs?
    if (any(is.na(edges)))
      stop("Some elements of -edges- were not able to be coerced as integers.")
    
    # Checking range
    ran <- range(edges)
    n   <- ran[2] + 1
    if (ran[1]!=0)
      stop("The minimum node id in the tree should be 0, it is ", ran[1], ".")
    
  }

  # Labels --------------------------------------------------------------------
  
  # Checking length of labels
  if (check.labels && (length(labels) != n)) {
    
    if (!is.vector(labels))
      stop("-labels- should be a vector.")
    
    if (!is.character(labels)) {
      warning("-labels-is not of class 'character'. Will be coerced.")
      labels <- as.character(labels)
    }
      
    stop("The length of -labels- should be equal to ", n, ". It is ",
         length(labels), ".")
    
  }
    
  
  # Offsprings -----------------------------------------------------------------
  
  # If no offspring was provided, then it is computed, and so no check should be
  # done
  
  if (!length(offspring)) {
    
    offspring       <- list_offspring(edges)
    check.offspring <- FALSE
    
  }
  
  if (check.offspring) {
    
    # Type list
    if (!is.list(offspring))
      stop("-offspring- should be a list.")
    
    # Length
    if (length(offspring) != n)
      stop("The length of -offspring- should be equal to ", n,". It is ",
           length(offspring), ".")
    
    # Subelements should be a list
    if (any(!sapply(offspring, is.vector)))
      stop("The elements of -offspring- should be vectors.")
    
    # Coercing into integer
    offspring <- lapply(offspring, as.integer)
    
    idx <- unlist(offspring, recursive = TRUE)
    
    # All should be integer
    if (any(is.na(idx)))
      stop("Some elements of -offspring- faild to be coerced to integer.")
    
    # All should be in the same range of ran
    ran_o <- range(idx)
    if (ran[1] > ran_o[1] | ran[2] < ran_o[2])
      stop("Some elements of -offspring- are pointing to elements out of range.")
    
  }
  
  # Branch lengths -------------------------------------------------------------
  if (!length(edge.length)) {
    
    edge.length       <- rep(1.0, nrow(edges))
    check.edge.length <- FALSE
    
  }
  
  if (check.edge.length) {
    
    if (!is.vector(edge.length))
      stop("-edge.length- should be a vector.")
    
    if (length(edge.length) != nrow(edges))
      stop("-edge.length- should have the same length as number of edges.")
    
    if (!is.numeric(edge.length))
      edge.length <- as.numeric(edge.length)
    
    if (any(!is.finite(edge.length))) 
      stop("Some -edge.length- show non-finite values.")
    
  }
  
  # Removing any other attribute that edges has --------------------------------
  
  mostattributes(edges) <- list(dim=dim(edges))
  attributes(edge.length)  <- NULL
  attributes(labels)    <- NULL
  attributes(offspring) <- NULL
  
  # Creating the final output
  structure(
    edges,
    class       = "po_tree",
    Nnode       = sum(sapply(offspring, length) > 0),
    edge.length = edge.length,
    labels      = labels,
    offspring   = offspring
  )
  
}

recode_as_po <- function(edges) {
  ans <- .recode_as_po(edges)
  new_po_tree(
    edges        = ans,
    labels       = attr(ans, "labels"), 
    check.edges  = FALSE,
    check.labels = FALSE
    )
}

#' @export
#' @param edge.length A numeric vector of length \code{nrow(edges)}. Branch
#' lengths.
#' @rdname as_po_tree
as_po_tree.default <- function(edges, edge.length=NULL, ...) {
  
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
  if (length(labels)) 
    labels <- labels[as.integer(attr(edges, "labels"))]
  
  new_po_tree(
    edges             = edges,
    labels            = labels,
    edge.length       = edge.length,
    offspring         = NULL, 
    check.edges       = FALSE,
    check.labels      = FALSE,
    check.edge.length = TRUE
  )

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
#' 
#' @details Plotting is done via \code{\link[ggtree:ggtree]{ggtree}} 
#' from the \pkg{ggtree} package (Bioconductor).
#' 
#' @return A list of class \code{aphylo} with the following elements:
#' \item{annotations}{A matrix of size \eqn{n\times P}{n * P} with leaf functional
#' annotations. These can be either 0, 1, or 9 indicating no function, function, and
#' no information respectively.}
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
    edges       = E
  )
  
}

as_aphylo <- function(
  annotations,
  edges,
  checks     = TRUE
  ) {
  
  if (checks) {
    
    # Checking class
    stopifnot(is.matrix(annotations))
    stopifnot(is.po_tree(edges))
    
    # Checking dimmensions
    d <- dim(annotations)
    stopifnot(length(unique(as.vector(edges))) == d[1])
    
  }
  
  structure(
    list(
      annotations = annotations,
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
  
  # Finding leafs
  is_leaf <- isleaf(E)
  
  # Figuring out labels
  tiplabs <- attr(x[["edges"]], "labels")[match(
    attr(E, "labels")[which(is_leaf)],
    attr(x[["edges"]], "labels")
    )]
  
  # Internal nodes labels
  nodelabs <- attr(x[["edges"]], "labels")[match(
    attr(E, "labels")[which(!is_leaf)],
    attr(x[["edges"]], "labels")
  )]
  
  attributes(E) <- list(dim=dim(E))
  
  structure(list(
    edge        = E,
    edge.length = rep(1, nrow(E)),
    tip.label   = as.character(unname(tiplabs)),
    Nnode       = sum(!is_leaf),
    node.label  = as.character(unname(nodelabs))
  ), class = "phylo")
  
}

#' @rdname as.apephylo
#' @export
as.apephylo.po_tree <- function(x, ...) {
  # Cleaning attributes
  attr(x, "offspring")  <- NULL
  
  # Recoding edgelist
  E <- as_ape_tree(x)
  
  is_leaf <- isleaf(E, from0 = FALSE)
  
  # Figuring out labels
  ids <- 0L:(length(attr(x, "labels")) - 1L)
  tiplabs <- attr(x, "labels")[match(attr(E, "labels")[which(is_leaf)], ids)]
  nodelabs <- attr(x, "labels")[match(attr(E, "labels")[which(!is_leaf)], ids)]
  
  attributes(E) <- list(dim=dim(E))
  
  structure(list(
    edge        = E,
    edge.length = attr(x, "edge.length"),
    tip.label   = as.character(unname(tiplabs)),
    node.label  = as.character(unname(nodelabs)),
    Nnode       = sum(!is_leaf)
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



#' Extract leaf labels 
#' @param x A phylogenetic tree.
#' @param ... Ignored.
#' @return \code{leafs} returns a character vector with the names of the leafs
#' @examples 
#' set.seed(1)
#' ans <- sim_tree(10)
#' leafs(ans)
#' @export
leafs <- function(x, ...) UseMethod("leafs")

#' @rdname leafs
#' @export
leafs.phylo <- function(x, ...) {
  x$tip.label
}

#' @rdname leafs
#' @export
leafs.po_tree <- function(x, ...) {
  labs <- attr(x, "labels")
  d <- fast_table_using_labels(x[,1], as.integer(names(labs)))
  
  unname(labs[which(d == 0)])
}

#' @rdname leafs
#' @export
leafs.aphylo <- function(x, ...) {
  leafs.po_tree(x$edges)
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
  
  ans <- lapply(1:ncol(object$annotations), function(i) table(object$annotations[,i]))
  
  if (!inherits(ans, "list"))
    ans <- list(ans[,1,drop=TRUE])
  
  ans <- do.call(
    rbind, 
    lapply(ans, function(x)
    data.frame(
      `0`  = unname(x["0"]),
      `1`  = unname(x["1"]),
      `NA` = unname(x["9"]),
      check.names = FALSE
      )
    )
  )
  ans[is.na(ans)] <- 0
  rownames(ans) <- colnames(object$annotations)
  cat("\nDistribution of functions:\n")
  print(ans)
  
  invisible(ans)
  
}

# This list sets the default plotting parameters when calling
# the plot.phylo function.
default.plot.phylo.params <- list(
  show.node.label = TRUE,
  show.tip.label  = TRUE,
  root.edge       = TRUE
)

# This is the actual function that does all the job setting the defaults.
set.default.plot.phylo.params <- function(dots) {
  env <- parent.frame()
  
  for (p in names(default.plot.phylo.params)) 
    if (!length(env[[dots]][[p]]))
      env[[dots]][[p]] <- default.plot.phylo.params[[p]]
}

#' @param x An object of class \code{aphylo}.
#' @param y Ignored.
#' @param ... Further arguments passed to the method.
#' @rdname as_po_tree
#' @export
plot.po_tree <- function(
  x, y=NULL, ...) {

  dots <- list(...)
  set.default.plot.phylo.params("dots")
  do.call(ape::plot.phylo, c(list(as.apephylo(x)), dots))
  
}

#' Relabel a tree (edgelist) using \pkg{ape}'s convention
#' 
#' This function takes an edgelist and recodes (relabels) the nodes following
#' \CRANpkg{ape}'s coding convention. \code{as_ape_tree} is the powerhorse of
#' \code{\link{as.apephylo}}.
#' 
#' @template parameters
#' @templateVar edges 1
#' @return An integer matrix of the same dimmension as \code{edges} with the following
#' aditional attribute:
#' \item{labels}{Named integer vector of size \code{n}. Original labels of the edgelist
#' where the first \code{n} are leaf nodes, \code{n+1} is the root node, and the reminder
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
  nodes <- c(leafs, roots, inner)
  
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
    class = c("matrix", "ape_tree")
  )
  
}


