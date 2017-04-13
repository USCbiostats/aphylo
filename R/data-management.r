#' Annotated Phylogenetic Tree
#' 
#' 
#' 
#' @param data_exper A data.frame with the experimental data
#' @param leafidvar A character scalar with the name of the leaf id variable
#' in \code{data_exper}.
#' @template parameters
#' @templateVar edges 1
#' @param nodeidvar A character scalar with the name of the node id variable
#' in \code{data_tree}.
#' @param parentidvar A character scalar with the name of the parent id variable
#' in \code{data_tree}.
#' @param funvars A character vector with the names of the function indicator
#' (0-1) variables in \code{data_exper}. If not provided, it is assumed that
#' all but \code{leafidvar} are function variables.
#' 
#' @details Plotting is done via \code{\link[ape:plot.phylo]{plot.phylo}} 
#' from the \CRANpkg{ape} package.
#' 
#' \code{as.phylo} requires loading the \code{ape} package.
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
  aphylo <- as_po_tree(edges)
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
#' \code{phylo} objects have several methods such as \code{plot}, \code{print}
#' and \code{summary}.
#'
#' @param x An object of class \code{phylo_tree} or \code{aphylo}.
#' @param ... Ignored.
#' @return An object of class \code{\link[ape:as.phylo]{phylo}}
#' @family Data management functions
#' @export
as.phylo <- function(x, ...) UseMethod("as.phylo")

#' @export
#' @rdname as.phylo
#' @method as.phylo default
as.phylo.default <- ape::as.phylo

#' @rdname new_aphylo
#' @method as.phylo aphylo
#' @export
as.phylo.aphylo <- function(x, ...) {
  
  # Recoding edgelist
  recoded_tree       <- with(x, recode_vertices(edges[,1], edges[,2]))
  
  structure(list(
    edge        = recoded_tree$edges,
    edge.length = rep(1, nrow(recoded_tree$edges)),
    tip.label   = sprintf("leaf%03i", 1:sum(recoded_tree$isleaf)),
    Nnode       = sum(!recoded_tree$isleaf)
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
  plot(as.phylo.aphylo(x), tip.color=tip.color, ...)
  
}

#' @method as.phylo phylo_tree
#' @rdname sim_tree
#' @export
as.phylo.phylo_tree <- function(x, ...) {
  
  # Recoding edgelist
  recoded_tree       <- recode_vertices(x[,1], x[,2])
  recoded_tree$edges <- recoded_tree$edges[,2:1]
  
  structure(list(
    edge        = recoded_tree$edges,
    edge.length = rep(1, nrow(recoded_tree$edges)),
    tip.label   = sprintf("leaf%03i", 1:sum(recoded_tree$isleaf)),
    Nnode       = sum(!recoded_tree$isleaf)
  ), class = "phylo")
  
}

#' @param x An object of class \code{phylo_tree}.
#' @param y Ignored.
#' @param tip.color See \code{\link[ape:plot.phylo]{plot.phylo}}
#' @param ... Further arguments passed to the method.
#' @rdname sim_tree
#' @export
plot.phylo_tree <- function(
  x, y=NULL, tip.color=NULL, ...) {
  # 
  # if (!length(tip.color)) {
  #   tip.color <- colors(9)[with(x, experiment[leaf_node,1,drop=TRUE])]
  # } 
  # 
  plot(as.phylo(x), tip.color=tip.color, ...)
  
}

#' Relabel a tree using \pkg{ape}'s convention
#' @param offspring Vector. Labels of offspring.
#' @param parent Vector. Labels of parents.
#' @param root_label Character scalar. Label of the root node.
#' @return A list with three components:
#' \item{edges}{Integer Matrix of size \code{length(offspring)*2}. Recoded edgelist}
#' \item{nodes}{Named character scalar of size \code{n}. Each node's id and label.}
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
#' ans0 <- sim_tree(30)
#' 
#' # In the ape package, nodes are labeled such that leafs are from 1 to n
#' # root node is n+1 and interior nodes can have whatever label.
#' ans1 <- recode_vertices(ans0[,1], ans0[,2])
#' 
#' head(ans0)
#' head(ans1)
#' 
#' @export
#' @family Data management functions
recode_vertices <- function(offspring, parent, root_label = "0") {
  
  # Tagging leaf (tip) nodes. These have degree 1
  edges   <- unname(cbind(offspring, parent))
  isleaf  <- table(edges) == 1
  
  nodes        <- sort(unique(as.vector(edges)))
  names(nodes) <- nodes
  isleaf       <- isleaf[match(names(isleaf), names(nodes))]
  
  # Leaf nodes must have ids 1:sum(isleaf), rootnode must have id sum(isleaf) + 1
  nodes[isleaf]                      <- 1:sum(isleaf)
  nodes["0"]                         <- sum(isleaf) + 1
  nodes[names(nodes) != 0 & !isleaf] <- (sum(isleaf) + 2):(length(isleaf))
  
  # Replacing elements in edges
  edges[] <- nodes[match(edges[], names(nodes))]
  edges   <- apply(edges, 2, as.integer)
  
  list(
    edges  = edges,
    nodes  = nodes,
    isleaf = isleaf
  )
  
}
