#' List each nodes' offspring
#' 
#' For each node in a tree, lists all its offspring.
#' 
#' @param x An object of class `phylo` or `aphylo`.
#' @return List of length `n` (total number of nodes).
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

#' @rdname list_offspring
#' @export
list_offspring.aphylo <- function(x) {
  x$offspring
}

#' @rdname list_offspring
#' @export
list_offspring.phylo <- function(x) {
  .list_offspring(x$edge, x$Nnode + length(x$tip.label))
}

map_ids_to_positions.aphylo_estimates <- function(ids_name, dat_name) {
  
  # Retrieving information about the tree 
  env    <- parent.frame()
  labels <- with(env[[dat_name]][["dat"]][["tree"]], c(tip.label, node.label))
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
      leafs    = 0L:(length(env[[dat_name]][["dat"]][["tree"]]$tip.label) - 1L),
      missings = {
        
        tmp <- with(env[[dat_name]][["dat"]], rbind(tip.annotation, node.annotation))
        which(!stats::complete.cases(tmp)) - 1L
        
      }
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


#' Annotated Phylogenetic Tree
#' 
#' The `aphylo` class tree holds both the tree structure represented as a
#' partially ordered phylogenetic tree, and node annotations. While annotations
#' are included for both leafs and inner nodes, the algorithms included in this
#' package only uses the leaf annotations.
#' 
#' @template parameters
#' @templateVar .tree 1
#' @templateVar .tip.annotation 1
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
#' ans <- new_aphylo(fakeexperiment[,2:3], faketree)
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
#' @export
new_aphylo <- function(
  tip.annotation,
  tree, 
  node.annotation = NULL
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
  
  offspring <- list_offspring(tree)
  
  structure(
    c(
      list(tree            = tree),
      list(tip.annotation  = tip.annotation),
      list(node.annotation = node.annotation),
      list(offspring       = offspring),
      list(pseq            = pseq),
      list(reduced_pseq    = reduce_pseq(pseq, rbind(tip.annotation, node.annotation), offspring))
      
    ),
    class = c("aphylo")
  )
}


#' Set of colors from dput(RColorBrewer::brewer.pal(7, "RdBu"), file = "")
#' @noRd
.aphyloColors <- c("#B2182B", "#EF8A62", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#67A9CF", 
                   "#2166AC")

#' Plot and print methods for `aphylo` objects
#' 
#' @param x An object of class `aphylo`.
#' @param y Ignored.
#' @param ... Further arguments passed to [ape::plot.phylo].
#' @param prop Numeric scalar between 0 and 1. Proportion of the device that the
#' annotations use in `plot.aphylo`.
#' @param rect.args List of arguments passed to [graphics::rect].
#' @name aphylo-methods
#' @details The `plot.aphylo` function is a wrapper of [ape::plot.phylo].
#' 
#' @export
#' @return In the case of `plot.aphylo`, `NULL`.
#' @family aphylo methods
#' @export
plot.aphylo <- function(
  x, y = NULL, prop = .15, 
  rect.args = list(), ...) {
  
  # Coercing into phylo
  phylo <- as.phylo(x)
  dots  <- list(...)
  
  # Some defaults
  if (!length(dots$cex))
    dots$cex <- .75
  if (!length(dots$show.node.label))
    dots$show.node.label <- TRUE
  if (!length(dots$font))
    dots$font <- 1
  if (!length(dots$main))
    dots$main <- "Annotated Phylogenetic Tree"
  
  if (length(dots$type) && dots$type != "phylogram")
    stop("Only `phylograph` is currently supported.")
  
  # Size of the device
  dev_size <- graphics::par("din")
  
  # How much space for the annotations
  labwidth <- dev_size[1]*prop

  op <- graphics::par(
    mai = graphics::par("mai")*c(0, 1, 1, 0) + c(labwidth, 0, 0, labwidth)
    )
  
  on.exit(graphics::par(op))
  do.call(graphics::plot, c(list(x=phylo), dots))
  
  # Capturing the parameters from the `ape` package
  plot_pars <- utils::getFromNamespace(".PlotPhyloEnv", "ape")
  
  tips <- with(plot_pars$last_plot.phylo, cbind(xx, yy))
  tips <- tips[1L:ape::Ntip(phylo),,drop=FALSE]
  
  yspacing <- range(tips[,2])
  yspacing <- (yspacing[2] - yspacing[1])/(nrow(tips) - 1)/2
  
  op2 <- graphics::par(
    mai = graphics::par("mai")*c(1,0,1,0) +
      c(0, dev_size[1]*(1-prop), 0, dev_size[1]*.025)
    )
  on.exit(graphics::par(op2), add=TRUE)
  
  graphics::plot.window(c(0, 1), range(tips[,2]), new=FALSE, xaxs = "i")
  
  nfun      <- ncol(x$tip.annotation)
  yran      <- range(tips[,2])
  
  graphics::rect(
    xleft   = -.1,
    ybottom = yran[1] - .1*yinch() - yspacing,
    xright  = 1.1,
    ytop    = yran[2] + .1*yinch() + yspacing,
    xpd     = NA,
    col     = "lightgray",
    border  = "gray"
  )
  
  for (f in 1:nfun) {

    rect.args$xleft   <- (f - 1)/nfun
    rect.args$ybottom <- tips[,2] - yspacing
    rect.args$xright  <- f/nfun
    rect.args$ytop    <- tips[,2] + yspacing
    rect.args$xpd     <- NA
    
    if (!length(rect.args$xpd)) rect.args$xpd <- NA
    if (!length(rect.args$col)) rect.args$col <- blue(x$tip.annotation[,f])
    if (!length(rect.args$border)) rect.args$border <- blue(x$tip.annotation[,f])
    if (!length(rect.args$lwd)) rect.args$lwd<-.5
    if (!length(rect.args$density)) rect.args$density <- 
      ifelse(x$tip.annotation[,f] == 9L, 10, NA)
    
    # Drawing rectangles
    do.call(graphics::rect, rect.args)
    
    # Adding function label
    graphics::text(
      x = (2*f - 1)/nfun/2 - 1/nfun/2,
      y = yran[1] - graphics::strheight(colnames(x$tip.annotation)[f], srt=45)*1.5 - yspacing,
      label = colnames(x$tip.annotation)[f],
      pos = 1,
      srt = 45,
      xpd = NA
    )
    
  }
  
  # Drawing a legend
  graphics::par(op2)
  graphics::par(mai = c(0,0,dev_size[2] - op2$mai[1], dev_size[1]*prop))
  graphics::plot.window(c(0,1), c(0,1))
  graphics::legend(
    "center",
    legend  = c("No function", "Function", "no information"),
    fill    = blue(c(0,1,9)),
    bty     = "n",
    density = c(NA, NA, 10),
    horiz   = TRUE
    )
  
  invisible(NULL)
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
APHYLO_DEFAULT_PLOT_PARAMS <- list(
  show.node.label = TRUE,
  show.tip.label  = TRUE,
  root.edge       = TRUE
)

# This is the actual function that does all the job setting the defaults.
set.default.plot.phylo.params <- function(dots) {
  dots <- as.character(match.call()$dots)
  env <- parent.frame()
  
  for (p in names(APHYLO_DEFAULT_PLOT_PARAMS)) 
    if (!length(env[[dots]][[p]]))
      env[[dots]][[p]] <- APHYLO_DEFAULT_PLOT_PARAMS[[p]]
}

