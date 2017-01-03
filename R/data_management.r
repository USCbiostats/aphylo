#' Find the offsprings of each node.
#' 
#' @param data_exper A data.frame with the experimental data
#' @param leafidvar A character scalar with the name of the leaf id variable
#' in \code{data_exper}.
#' @param data_tree A data.frame with the tree.
#' @param nodeidvar A character scalar with the name of the node id variable
#' in \code{data_tree}.
#' @param parentidvar A character scalar with the name of the parent id variable
#' in \code{data_tree}.
#' 
#' @details In the case of the plot method, by defult layout is set to
#' be \code{\link[igraph:layout_with_sugiyama]{layout_with_sugiyama}}
#' with parameter \code{maxiter=200}.
#' 
#' @return A list of length \eqn{n} with relative position of offsprings
#' of each node, with respect to \code{data_exper}, starting from 0.
#' @examples 
#' # Loading data
#' data(experiment)
#' data(tree)
#' 
#' ans <- get_offsprings(
#'     experiment, "LeafId", 
#'     tree, "NodeId", "ParentId"
#' )
#' 
#' str(ans)
#' 
#' # We can visualize it
#' plot(ans, vertex.size=5, vertex.label=NA)
#' @export
get_offsprings <- function(
  data_exper, 
  leafidvar,
  data_tree,
  nodeidvar,
  parentidvar
  ) {
  
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
      offsprings  = ans,
      noffsprings = sapply(ans, length),
      edgelist    = as.matrix(data_tree[,c(nodeidvar, parentidvar), drop=FALSE])
    ),
    class = "phylo_offsprings"
  )
}
