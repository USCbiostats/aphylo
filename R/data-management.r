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
#' @details In the case of the plot method, by defult layout is set to
#' be \code{\link[igraph:layout_with_sugiyama]{layout_with_sugiyama}}
#' with parameter \code{maxiter=200}.
#' 
#' @return A list of length \eqn{n} with relative position of offspring
#' of each node, with respect to \code{data_exper}, starting from 0.
#' @examples 
#' # Loading data
#' data(experiment)
#' data(tree)
#' 
#' ans <- get_offspring(
#'     experiment, "LeafId", 
#'     tree, "NodeId", "ParentId"
#' )
#' 
#' str(ans)
#' 
#' # We can visualize it
#' plot(ans, vertex.size=5, vertex.label=NA)
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
  
  # Checking who is parent
  structure(
    list(
      experiment = unname(as.matrix(data_exper[,funvars])),
      fun_names  = funvars, 
      added      = data_exper[["added"]],
      offspring  = ans,
      noffspring = sapply(ans, length),
      edgelist   = as.matrix(data_tree[,c(nodeidvar, parentidvar), drop=FALSE])
    ),
    class = "phylo_offspring"
  )
}

#' Convert \code{phylo_offspring} to \code{igraph} objects.
#' @param x An object of class \code{phylo_offspring}
#' @export
phylo_offspring_to_igraph <- function(x) {
  
  # Checking class
  if (!inherits(x, "phylo_offspring"))
    stop("-x- must be of class -phylo_offspring-.")
  
  # Retrieving the edgelist
  edgelist   <- x[["edgelist"]]
  edgelist[] <- as.character(edgelist)
  
  # Converting into igraph object
  igraph::graph_from_edgelist(edgelist)
}

#' @rdname phylo_offspring_to_igraph
#' @export
phylo_o_to_phylo <- function(x) {
  
  # Checking class
  if (!inherits(x, "phylo_offspring"))
    stop("-x- must be of class -phylo_offspring-.")
  
  # Recoding edgelist
  graph <- phylo_offspring_to_igraph(x)
  graph <- data.frame(igraph::as_edgelist(graph), stringsAsFactors=TRUE)
  
  
  
  graph
  
}

