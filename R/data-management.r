#' List each nodes' offspring
#' 
#' For each node in a tree, lists all its offspring.
#' 
#' @param x Either a matrix or an object of class `phylo`
#' @return A named list of length `n` (total number of nodes). The names of the
#' list correspond to the ids of the nodes in `x`.
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

list_offspring.aphylo <- function(x) {
  x$offspring
}

list_offspring.phylo <- function(x) {
  ids <- 1L:(x$Nnode + length(x$tip.label))
  
  structure(
    .list_offspring(x$edge, length(ids)),
    names = ids
  )
  
}

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

map_ids_to_positions.aphylo_estimates <- function(ids_name, dat_name) {
  
  # Retrieving information about the tree 
  env <- parent.frame()
  labels <- attr(env[[dat_name]][["dat"]][["edges"]], "labels")
  n      <- length(labels)
  
  # Mapping the ids to the positions -------------------------------------------
  if (is.numeric(env[[ids_name]])) {
    
    # Coercing into integer
    env[[ids_name]] <- as.integer(env[[ids_name]]) 
    
  } else if (length(env[[ids_name]]) == 1 && (env[[ids_name]] %in% c("all", "leafs", "missings"))) {
    
    # Matching the description
    env[[ids_name]] <- switch(
      env[[ids_name]], 
      all      = 0L:(n - 1L),
      leafs    = which(isleaf(env[[dat_name]][["dat"]][["edges"]])) - 1L,
      missings = which(
        isleaf(env[[dat_name]][["dat"]][["edges"]]) &
          apply(env[[dat_name]][["dat"]][["annotations"]], 1, function(a) any(a == 9)*1L) > 0
      ) - 1L
    )
    
  } else if (is.character(env[[ids_name]])) {
    
    # Matching the labels
    env[[ids_name]] <- as.integer(match(env[[ids_name]], labels)) - 1L
    
  } else
    stop("Unsupported type of ids: ", paste0(env[[ids_name]], collapse=", "))
  
  # Checking everything is in order --------------------------------------------
  
  # All complete
  if (any(!is.finite(env[[ids_name]])))
    stop("Some ids are either NAs or Inf.")
  
  if (any(env[[ids_name]] >= n))
    stop("Out of range: Some ids are greater than (n-1). Ids go from 0 to n-1.")
  
  if (any(env[[ids_name]] < 0L))
    stop("Out of range: Some ids are greater less than 0. Ids go from 0 to n-1.")
  
}

map_ids_to_positions.po_tree <- function(ids_name, edges_name) {
  
  # Retrieving information about the tree 
  env <- parent.frame()
  labels <- attr(env[[edges_name]], "labels")
  n      <- length(labels)
  
  # Sorting and unique
  env[[ids_name]] <- sort(unique(env[[ids_name]]))
  
  # Mapping the ids to the positions -------------------------------------------
  if (is.numeric(env[[ids_name]])) {
    
    # Coercing into integer
    env[[ids_name]] <- as.integer(env[[ids_name]]) 
    
  } else if (length(env[[ids_name]]) == 1 && (env[[ids_name]] %in% c("all", "leafs"))) {
    
    # Matching the description
    env[[ids_name]] <- switch(
      env[[ids_name]], 
      all      = 0L:(n - 1L),
      leafs    = which(isleaf(env[[edges_name]])) - 1L
    )
    
  } else if (is.character(env[[ids_name]])) {
    
    # Matching the labels
    env[[ids_name]] <- as.integer(match(env[[ids_name]], labels)) - 1L
    
  } else
    stop("Unsupported type of -ids-: ", paste0(env[[ids_name]], collapse=", "))
  
  # Checking everything is in order --------------------------------------------
  
  # All complete
  if (any(!is.finite(env[[ids_name]])))
    stop("Some ids are either NAs or Inf.")
  
  if (any(env[[ids_name]] >= n))
    stop("Out of range: Some ids are greater than (n-1). Ids go from 0 to n-1.")
  
  if (any(env[[ids_name]] < 0L))
    stop("Out of range: Some ids are greater less than 0. Ids go from 0 to n-1.")
  
}

#' Recodes an edgelist as a Partially Ordered Tree
#' 
#' The function [new_aphylo()] uses this function to make sure that
#' the edgelist that is passed makes a partial order. This is a requirement
#' for the peeling algorithm, which is used explicitly in the `LogLike`
#' function.
#' 
#' @details
#' The recoded edgelist is such that in all rows the first element, parent
#' node, as a label that is less than the second element, the offspring, a
#' partial order.
#' 
#' @return A matrix of the same dimension as `edges`, an edgelist, recoded
#' to form a partial order. Besides of been of class `matrix`, the resulting
#' object is also of class `po_tree` and has an aditional attributes:
#' \item{Nnode}{Integer scalar. Number of leaf nodes.}
#' \item{edge.length}{Numeric vector of length `nrow(edges)`. Branch lengths.}
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
#' as.phylo(potree)
#'
#' @aliases po_tree 
#' @export
as_po_tree <- function(edges, ...) UseMethod("as_po_tree")

#' @export
#' @rdname as_po_tree
#' @details `as.po_tree` is an alias for `as_po_tree`
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
#' @return In the case of `is.po_tree`, `TRUE` if is an object of
#' class `po_tree`, `FALSE` otherwise.
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
  edge.length     = NULL,
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
    class       = c("po_tree", "matrix"),
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
#' @param edge.length A numeric vector of length `nrow(edges)`. Branch
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
#' The `aphylo` class tree holds both the tree structure represented as a
#' partially ordered phylogenetic tree, and node annotations. While annotations
#' are included for both leafs and inner nodes, the algorithms included in this
#' package only uses the leaf annotations.
#' 
#' @template parameters
#' @templateVar edges 1
#' @templateVar annotationslab 1
#' 
#' @details Plotting is done via [ggtree:ggtree::ggtree()] 
#' from the \pkg{ggtree} package (Bioconductor).
#' 
#' @return A list of class `aphylo` with the following elements:
#' \item{tree}{An object of class [phylo][ape::read.tree].}
#' \item{tip.annotation}{An integer matrix. Tip (leaf) nodes annotations.}
#' \item{node.annotation}{An integer matrix (optional). Internal nodes
#'   annotations.}
#' \item{offspring}{A list. List of offspring of each node.}
#' \item{pseq}{Integer vector. The pruning sequence (postorder).}
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


#' Checks whether the annotations are in the correct format.
#' @param x Some node level annotations
#' @details 
#' The criteria is:
#' - Either a matrix or a data frame
#' - Checks the length
#' - Should have at least 2 columns
#' - No duplicates
#' - Should be numeric
#' - 0, 1, 9, or NA.
#' - NAs are replaces with 9s
#' - Set colnames
#' @noRd
check_annotations <- function(x) {

  # Nothing to check if it is NULL
  if (!length(x))
    return(x)
  
  # If it is a vector, then coerce it into a one column matrix
  if (is.vector(x) && !is.list(x))
    x <- cbind(x)
  
  
  # Checking class
  if (!inherits(x, c("data.frame", "matrix")))
    stop("-annotations- must be either an object of class -data.frame- or -matrix- (is of class ",
         class(x),").")

  test <- apply(x, 2, typeof)
  test <- which(!sapply(test, `%in%`, table = c("integer", "double")))
  if (length(test))
    stop("All annotation columns (but the first one) must be either -integer- or -double-.",
         " The following columns are not: ", paste(test, collapse=", "), ".")
  
  node_annotations <- apply(x, 2, as.integer)
  
  # Checking values
  test <- which(apply(node_annotations, 1, function(y) any(!(y %in% c(0, 1, 9, NA)))))
  if (length(test))
    stop("The following rows of -annotations- have values different from c(0, 1, 9, NA): ",
         paste(test, collapse=", "), ".")
  
  # Replacing NAs
  node_annotations[is.na(node_annotations)] <- 9L
  
  # Checking colnames
  fun_names <- colnames(node_annotations)
  
  if (!length(fun_names))
    fun_names <- sprintf("fun%03i", P)
  
  # Setting dimnames
  dimnames(node_annotations) <- list(
    1L:nrow(node_annotations),
    fun_names
  )
  
  # Actual annotations
  node_annotations
  
}

#' @rdname aphylo-class
#' @export
new_aphylo <- function(
  tip.annotation,
  node.annotation = NULL,
  tree
  ) {
  
  # Coercing tree to a phylo object
  tree <- ape::as.phylo(tree)
  
  # Checking annotations
  tip.annotation  <- check_annotations(tip.annotation)
  node.annotation <- check_annotations(node.annotation)
  
  # Are all annotations in the edgelist
  if (nrow(tip.annotation) != length(tree$tip.label))
    stop("The number of `tip.annotation` differs with the number of tips in `tree`.",
         call. = FALSE)
  
  if (length(node.annotation) && (nrow(node.annotation) != tree$Nnode))
    stop("The number of `node.annotation` differs with the number of internal nodes in `tree`.",
         call. = FALSE)
  
  # Returning
  as_aphylo(
    tip.annotation  = tip.annotation,
    node.annotation = node.annotation,
    tree             = tree,
    checks           = FALSE
  )
  
}

#' For internal use only: Creates the aphylo class object
#' @noRd
#' 
as_aphylo <- function(
  tip.annotation,
  node.annotation,
  tree,
  checks     = TRUE
  ) {
  
  if (checks) {
    
    # Checking class
    stopifnot(is.matrix(tip.annotation))
    stopifnot(inherits(tree, "phylo"))
    
    # Checking dimmensions
    stopifnot(nrow(tip.annotation) == length(tree$tip.label))
    if (length(node.annotation)) {
      stopifnot(is.matrix(node.annotation))
      stopifnot(nrow(node.annotation) == tree$Nnode)
    } 
      
  }
  
  # PATCH FIX: FOR NOW WE NEED TO HAVE SOMETHING FOR NODES!
  if (!length(node.annotation))
    node.annotation <- matrix(
      9L, nrow = tree$Nnode,
      ncol = ncol(tip.annotation),
      dimnames = list(
        (nrow(tip.annotation) + 1):(nrow(tip.annotation) + tree$Nnode),
        colnames(tip.annotation)
      )
    )
  
  # Calculating the pruning sequence
  pseq <- ape::postorder(tree)
  pseq <- c(tree$edge[pseq, 2], length(tree$tip.label) + 1L)
  
  structure(
    c(
      list(tree = tree),
      list(tip.annotation = tip.annotation),
      list(node.annotation = node.annotation),
      list(offspring = list_offspring(tree)),
      list(pseq      = pseq)
    ),
    class = c("aphylo")
  )
}


#' Plot and print methods for `aphylo` objects
#' 
#' @param x An object of class `aphylo`.
#' @param y Ignored.
#' @param geom.tiplab.args Further arguments passed to [ggtree:ggtree::ggtree()]
#' @param gheatmap.args Further arguments passed to [ggtree:ggtree::ggtree()]
#' @param scale.fill.args Further arguments passed to [ggtree:ggtree::ggtree()]
#' @param ... Further arguments passed to the method.
#' @name aphylo-methods
#' @details The `plot.aphylo` function is a wrapper of [ggtree:ggtree::ggtree()]
#' that creates a visualization as follows:
#' \enumerate{
#'   \item Retrieve the annotations from the [=aphylo-class::aphylo()] object
#'   \item Create a `ggtree` map adding [ggtree:geom_tiplab::geom_tiplab()]s
#'   \item Use [ggtree:gheatmap::gheatmap()] to add a heatmap.
#'   \item Set the colors using [ggplot2:scale_fill_manual::scale_fill_manual()]
#' }
#' 
#' @export
#' @return In the case of `plot.aphylo`, an object of class `c("ggtree", "gg", "ggplot")`
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
  if (!("colnames_angle" %in% names(gheatmap.args))) gheatmap.args$colnames_angle <- 45
  
  if (!("breaks" %in% names(scale.fill.args))) scale.fill.args$breaks <- c("0", "1", "9")
  if (!("values" %in% names(scale.fill.args))) scale.fill.args$values <- c("gray", "steelblue", "white")
  if (!("labels" %in% names(scale.fill.args))) scale.fill.args$labels <- c("No function", "Function", "N/A")
  
  # Coercing into a tree
  tree <- as.phylo(x)
  
  # Retrieving the annotations
  A   <- x$tip.annotation
  A[] <- as.character(A)
  
  # Matching positions
  A   <- A[match(tree$tip.label, rownames(A)),,drop=FALSE]
  
  # Creating the mapping
  p <- ggtree::ggtree(tree) +
    do.call(ggtree::geom_tiplab, geom.tiplab.args)
  
  # Adding the functions
  p <- do.call(ggtree::gheatmap, c(list(p=p, data=A),gheatmap.args)) +
    do.call(ggplot2::scale_fill_manual, scale.fill.args)
  
  p
}



#' Extract leaf labels 
#' @param x A phylogenetic tree.
#' @param ... Ignored.
#' @return `leafs` returns a character vector with the names of the leafs
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
#' @return In the case of `print.aphylo`, the same `aphylo` object
#' invisible.
print.aphylo <- function(x, ...) {
  
  print(x$tree)
  
  cat("\n Tip (leafs) annotations:\n")
  
  print(utils::head(x$tip.annotation))
  
  if (nrow(x$tip.annotation) > 6)
    cat("\n...(", nrow(x$tip.annotation) - 6, " obs. omitted)...\n\n", sep="")
  
  # Printing node annotations
  if (length(x$node.annotation)) {
    cat("\n Internal node annotations:\n")
    
    print(utils::head(x$node.annotation))
    
    if (nrow(x$node.annotation) > 6)
      cat("\n...(", nrow(x$node.annotation) - 6, " obs. omitted)...\n\n", sep="")
  } else {
    cat("\nNo annotations for internal nodes.")
  }
  

  invisible(x)
}

#' @export
#' @param object An object of class `aphylo`.
#' @rdname aphylo-methods
summary.aphylo <- function(object, ...) {
  
  for (x in c("tip.annotation", "node.annotation")) {
    ans <- lapply(1:ncol(object[[x]]), function(i) table(object[[x]][,i]))
    
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
    rownames(ans) <- colnames(object[[x]])
    cat("\nDistribution of functions in", x, ":\n")
    print(ans)
  }
  
  
  
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
  dots <- as.character(match.call()$dots)
  env <- parent.frame()
  
  for (p in names(default.plot.phylo.params)) 
    if (!length(env[[dots]][[p]]))
      env[[dots]][[p]] <- default.plot.phylo.params[[p]]
}

#' @param x An object of class `aphylo`.
#' @param y Ignored.
#' @param ... Further arguments passed to the method.
#' @rdname as_po_tree
#' @export
plot.po_tree <- function(
  x, y=NULL, ...) {

  dots <- list(...)
  set.default.plot.phylo.params(dots)
  do.call(ape::plot.phylo, c(list(as.phylo(x)), dots))
  
}


