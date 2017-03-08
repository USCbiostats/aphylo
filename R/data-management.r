#' Find the offspring of each node.
#' 
#' @param data_exper A data.frame with the experimental data
#' @param leafidvar A character scalar with the name of the leaf id variable
#' in \code{data_exper}.
#' @param data_tree A data.frame with the tree.
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
#' # Example with lots of data -------------------------------------------------
#' # Loading data
#' data(experiment)
#' data(tree)
#' 
#' ans <- get_offspring(
#'     experiment, "LeafId", 
#'     tree, "NodeId", "ParentId"
#' )
#' 
#' # We can visualize it
#' plot(ans)
#' 
#' # Example with less data ----------------------------------------------------
#' 
#' ans <- get_offspring(
#'     fakeexperiment, "LeafId", 
#'     faketree, "NodeId", "ParentId"
#' )
#'  
#' # We can visualize it
#' plot(ans)
#' @export
#' @family Data management functions
get_offspring <- function(
  data_exper,
  leafidvar,
  data_tree,
  nodeidvar,
  parentidvar,
  funvars=NULL
  ) {
  
  # Identifying function variables
  if (!length(funvars))
    funvars <- colnames(data_exper)[colnames(data_exper) != leafidvar]
  
  # Completing -data_exper- with the tree
  ids        <- unique(unlist(data_tree[,c(nodeidvar, parentidvar)]))
  test       <- which(!(ids %in% data_exper[[leafidvar]]))
  if (length(test)) {
    
    # Adding the tag variable (default FALSE)
    data_exper[["added"]] <- FALSE
    
    # Creating the rows to add
    newrows              <- data_exper[rep(1, length(test)),,drop=FALSE]
    newrows[,funvars]    <- 9
    newrows[[leafidvar]] <- ids[test]
    newrows[["added"]]   <- TRUE
    
    # Adding the rows
    data_exper <- rbind(data_exper, newrows)
    rownames(data_exper) <- NULL
  }
  
  # Checking Data sorting
  data_exper <- data_exper[
    order(data_exper[[leafidvar]], decreasing = FALSE),,drop=FALSE
    ]
  
  # Finding 
  ans <- lapply(data_exper[[leafidvar]], function(x) {
    x <- data_tree[[nodeidvar]][which(data_tree[[parentidvar]] == x)]
    x <- match(x, data_exper[[leafidvar]])
    x[!is.na(x)]
  })
  
  # Substracting one so we canuse it in C++
  ans <- lapply(ans, function(x) {
    if (length(x)) return(x-1)
    else return(x)
  })
  
  # Returning
  as_phylo_offspring(
    experiment = unname(as.matrix(data_exper[,funvars])),
    fun_names  = funvars, 
    added      = data_exper[["added"]],
    offspring  = ans,
    noffspring = sapply(ans, length),
    edgelist   = unname(as.matrix(data_tree[,c(nodeidvar, parentidvar), drop=FALSE]))
  )
  
  # structure(
  #   list(
  #     experiment = unname(as.matrix(data_exper[,funvars])),
  #     fun_names  = funvars, 
  #     added      = data_exper[["added"]],
  #     offspring  = ans,
  #     noffspring = sapply(ans, length),
  #     edgelist   = unname(as.matrix(data_tree[,c(nodeidvar, parentidvar), drop=FALSE]))
  #   ),
  #   class = "phylo_offspring"
  # )
}

as_phylo_offspring <- function(
  experiment,
  fun_names,
  added,
  offspring,
  noffspring,
  edgelist,
  checks = FALSE) {
  
  structure(
    list(
      experiment = experiment,
      fun_names  = fun_names, 
      added      = added,
      offspring  = offspring,
      noffspring = noffspring,
      edgelist   = edgelist
    ),
    class = "phylo_offspring"
  )
}

#' Coercing into \code{phylo} objects of the \pkg{ape} package.
#'
#' \code{phylo} objects have several methods such as \code{plot}, \code{print}
#' and \code{summary}.
#'
#' @param x An object of class \code{phylo_tree} or \code{phylo_offspring}.
#' @param ... Ignored.
#' @return An object of class \code{\link[ape:as.phylo]{phylo}}
#' @family Data management functions
#' @export
as.phylo <- function(x, ...) UseMethod("as.phylo")

#' @export
#' @rdname as.phylo
#' @method as.phylo default
as.phylo.default <- ape::as.phylo

#' @rdname get_offspring
#' @method as.phylo phylo_offspring
#' @export
as.phylo.phylo_offspring <- function(x, ...) {
  
  # Recoding edgelist
  recoded_tree       <- with(x, recode_vertices(edgelist[,1], edgelist[,2]))
  recoded_tree$edges <- recoded_tree$edges[,2:1]
  
  structure(list(
    edge        = recoded_tree$edges,
    edge.length = rep(1, nrow(recoded_tree$edges)),
    tip.label   = sprintf("leaf%03i", 1:sum(recoded_tree$isleaf)),
    Nnode       = sum(!recoded_tree$isleaf)
  ), class = "phylo")
  
}



#' Plot and print methods for \code{phylo_offspring} objects
#' @param x An object of class \code{phylo_offspring}.
#' @param y Ignored.
#' @param tip.color See \code{\link[ape:plot.phylo]{plot.phylo}}
#' @param ... Further arguments passed to the method.
#' @rdname get_offspring
#' @export
plot.phylo_offspring <- function(
  x, y=NULL, tip.color=NULL, ...) {
  # 
  # if (!length(tip.color)) {
  #   tip.color <- colors(9)[with(x, experiment[leaf_node,1,drop=TRUE])]
  # } 
  # 
  plot(as.phylo.phylo_offspring(x), tip.color=tip.color, ...)
  
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
