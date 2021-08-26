#' Annotated Phylogenetic Tree
#' 
#' The `aphylo` class tree holds both the tree structure represented as a
#' partially ordered phylogenetic tree, and node annotations. While annotations
#' are included for both leafs and inner nodes, the algorithms included in this
#' package only uses the leaf annotations.
#' 
#' @template parameters
#' @templateVar .tree 1
#' @templateVar .types 1
#' @templateVar .tip.annotation 1
#' 
#' @return A list of class `aphylo` with the following elements:
#' \item{tree}{An object of class [phylo][ape::read.tree].}
#' \item{tip.annotation}{An integer matrix. Tip (leaf) nodes annotations.}
#' \item{node.annotation}{An integer matrix (optional). Internal nodes
#'   annotations.}
#' \item{offspring}{A list. List of offspring of each node.}
#' \item{pseq}{Integer vector. The pruning sequence (postorder).}
#' \item{reduced_pseq}{Integer vector. The reduced version of `pseq`.}
#' \item{Ntips.annotated}{Integer. Number of tips with annotations.}
#' \item{tip.type}{Binary of length [Ntip()]. 0 means duplication and 1 speciation.}
#' \item{tip.type}{Binary of length [Nnode()]. 0 means duplication and 1 speciation.}
#' 
#' @examples 
#' # A simple example ----------------------------------------------------------
#' 
#' data(fakeexperiment)
#' data(faketree)
#' ans <- new_aphylo(fakeexperiment[,2:3], tree = as.phylo(faketree))
#'  
#' # We can visualize it
#' plot(ans)
#' @family Data management functions
#' @family aphylo methods
#' @name aphylo-class
#' @aliases aphylo
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
    fun_names <- sprintf("fun%03i", ncol(node_annotations))
  
  # Setting dimnames
  dimnames(node_annotations) <- list(
    1L:nrow(node_annotations),
    fun_names
  )
  
  # Actual annotations
  node_annotations
  
}

#' @rdname aphylo-class
#' @param ... Further argmuents passed to the method.
#' @export
new_aphylo <- function(tree, tip.annotation, ...) UseMethod("new_aphylo")

#' @rdname aphylo-class
#' @export
new_aphylo.phylo <- function(
  tree,
  tip.annotation,
  node.annotation = NULL,
  tip.type        = NULL,
  node.type       = NULL,
  ...
  ) {
  
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
  
  if (is.null(tip.type)) {
    tip.type <- integer(ape::Ntip(tree))
  } else {
    
    tip.type <- as.vector(tip.type)
    if (length(tip.type) != ape::Ntip(tree))
      stop(
        "The provided -tip.type- has not the same length as the number of tips",
        " in the tree. The actual length should be ", ape::Ntip(tree),
           call. = FALSE)
    
  }
  
  if (is.null(node.type)) {
    node.type <- integer(ape::Nnode(tree))
  } else {
    
    node.type <- as.vector(node.type)
    if (length(node.type) != ape::Nnode(tree))
      stop(
        "The provided -node.type- has not the same length as the number of tips",
        " in the tree. The actual length should be ", ape::Nnode(tree),
           call. = FALSE)
    
  }
  
  # Returning
  as_aphylo(
    tip.annotation  = tip.annotation,
    node.annotation = node.annotation,
    tree            = tree,
    tip.type        = tip.type,
    node.type       = node.type,
    checks          = FALSE
  )
  
}

#' For internal use only: Creates the aphylo class object
#' @noRd
#' 
as_aphylo <- function(
  tip.annotation,
  node.annotation,
  tree,
  tip.type,
  node.type,
  checks = TRUE
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
    
    # The nunmber of types should be of the same size of data points
    if (length(tip.type))
      stopifnot(length(tip.type) == ape::Ntip(tree))
    if (length(node.type))
      stopifnot(length(node.type) == ape::Nnode(tree))
    
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
  
  # If empty, then we fill it. Otherwise the likelihood function will fail (badly)
  if (!length(tip.type))
    tip.type <- integer(ape::Ntip(tree))
  if (!length(node.type))
    node.type <- integer(ape::Nnode(tree))
  
  # Listing offpring
  offspring <- list_offspring(tree)
  
  # Calculating the pruning sequence
  pseq         <- ape::postorder(tree)
  pseq         <- c(tree$edge[pseq, 2], length(tree$tip.label) + 1L)
  pseq_reduced <- reduce_pseq(pseq, rbind(tip.annotation, node.annotation), offspring)
  
  structure(
    c(
      list(tree            = tree),
      list(tip.annotation  = tip.annotation),
      list(node.annotation = node.annotation),
      list(offspring       = offspring),
      list(pseq            = pseq),
      list(reduced_pseq    = pseq_reduced),
      list(Ntips.annotated = length(intersect(1:nrow(tip.annotation), pseq_reduced))),
      list(tip.type        = tip.type),
      list(node.type       = node.type)
    ),
    class = c("aphylo")
  )
}


#' @export
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
summary.aphylo <- function(object, ...) {
  
  ans <- list()
  for (x in c("tip.annotation", "node.annotation")) {
    ans[[x]] <- lapply(1:ncol(object[[x]]), function(i) table(object[[x]][,i]))
    
    if (!inherits(ans[[x]], "list"))
      ans[[x]] <- list(ans[[x]][,1,drop=TRUE])
    
    ans[[x]] <- do.call(
      rbind, 
      lapply(ans[[x]], function(x)
        data.frame(
          `0`  = unname(x["0"]),
          `1`  = unname(x["1"]),
          `NA` = unname(x["9"]),
          check.names = FALSE
        )
      )
    )
    
    ans[[x]][is.na(ans[[x]])] <- 0
    rownames(ans[[x]]) <- colnames(object[[x]])
    cat("\nDistribution of functions in", x, ":\n")
    print(ans[[x]])
  }
  
  invisible(ans)
  
}


# # This list sets the default plotting parameters when calling
# # the plot.phylo function.
# APHYLO_DEFAULT_PLOT_PARAMS <- list(
#   show.node.label = TRUE,
#   show.tip.label  = TRUE,
#   root.edge       = TRUE
# )
# 
# # This is the actual function that does all the job setting the defaults.
# set.default.plot.phylo.params <- function(dots) {
#   dots <- as.character(match.call()$dots)
#   env <- parent.frame()
#   
#   for (p in names(APHYLO_DEFAULT_PLOT_PARAMS)) 
#     if (!length(env[[dots]][[p]]))
#       env[[dots]][[p]] <- APHYLO_DEFAULT_PLOT_PARAMS[[p]]
# }


is.aphylo <- function(x)      inherits(x, "aphylo")
is.multiAphylo <- function(x) inherits(x, "multiAphylo") 

