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
  
  # Tagging leaf (tip) nodes. These have degree 1
  edges  <- unname(as.matrix(data_tree[,c(nodeidvar, parentidvar), drop=FALSE]))
  isleaf <- table(edges) == 1
  
  nodes        <- sort(unique(as.vector(edges)))
  names(nodes) <- nodes
  isleaf       <- isleaf[match(names(isleaf), names(nodes))]
  
  # Leaf nodes must have ids 1:sum(isleaf), rootnode must have id sum(isleaf) + 1
  nodes[isleaf]                      <- 1:sum(isleaf)
  nodes["0"]                         <- sum(isleaf) + 1
  nodes[names(nodes) != 0 & !isleaf] <- (sum(isleaf) + 2):(length(isleaf))
  
  # Replacing elements in edges
  edges[] <- as.integer(nodes[match(edges[], names(nodes))])
  
  # Checking who is parent
  structure(
    list(
      experiment = unname(as.matrix(data_exper[,funvars])),
      fun_names  = funvars, 
      added      = data_exper[["added"]],
      offspring  = ans,
      noffspring = sapply(ans, length),
      edgelist   = edges,
      nodes      = nodes,
      leaf_node  = isleaf
    ),
    class = "phylo_offspring"
  )
}

#' @rdname get_offspring
#' @method as.phylo phylo_offspring
#' @export
as.phylo.phylo_offspring <- function(x, ...) {
  
  # Checking class
  if (!inherits(x, "phylo_offspring"))
    stop("-x- must be of class -phylo_offspring-.")
  
  # Recoding edgelist
  edges <- x$edgelist[,2:1]
  
  structure(list(
    edge        = unname(edges),
    edge.length = rep(1, nrow(edges)),
    tip.label   = sprintf("leaf%03i", 1:sum(x$leaf_node)),
    Nnode       = sum(!x$leaf_node)
  ), class = "phylo")
  
}

#' #' @export
#' #' @rdname get_offspring
#' as.phylo <- ape::as.phylo

# setGeneric("as.phylo")
# 
# O2 <- get_offspring(fakeexperiment, "LeafId", faketree, "NodeId", "ParentId")
# ans <- phylo_o_to_phylo(O2)
# collapse.singles(ans)
# plot(ans)
