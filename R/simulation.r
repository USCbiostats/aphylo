#' Random tree generation
#' 
#' By randomly choosing pairs of vertices, this function generates random
#' trees from bottom to top such that the labels of the nodes follow a 
#' partial order in their parent-offspring relation, with the parent
#' always having a lower idlabel than the offspring.
#' 
#' @param n Integer scalar. Number of leaf nodes.
#' @param edge.length A Function. Used to set the length of the edges.
#' 
#' @details The algorithm was implemented as follows
#' 
#' \enumerate{
#'   \item Initialize `left =[n*2 - 2,...,(n-1)]` and `m = n*2 - 2`, and
#'         initialize the vectors `parent` and `offspring` to be
#'         empty.
#'   \item While `length(left) > 1` do:
#'   \enumerate{
#'     \item Randomly choose a pair `(i, j)` from `left`
#'     \item Add `leaf(i)`, `leaf(j)` to the tail of `offspring`,
#'     \item Decrease `m` by 1, and add it two times to the tail of
#'     `parent`.
#'     \item Remove `(i,j)` from `leaf` and add `m` to its tail.
#'     \item next
#'   }
#'   
#' }
#' 
#' The [ape:rtree::rtree()] function in the \pkg{ape} package is similar,
#' although the big difference is in the way the labels are stablished. This later
#' point is crucial for both \pkg{phylogenetic} and \pkg{ape} as is a key feature
#' in some (most) of its routines.
#' 
#' @return An matrix of size `(n*2 - 2)*2`, an edgelist, with `n*2-1` nodes.
#' This, a Directed Acyclic Graph (DAG), as classes `matrix` and `po_tree`.
#' With the following additional attributes:
#' \item{offspring}{A list of size `n*2 - 1` listing node ith's offspring if any.}
#' 
#' @examples
#' # A very simple example ----------------------------------------------------
#' set.seed(1223)
#' newtree <- sim_tree(50)
#' 
#' plot(as.apephylo(newtree))
#' 
#' # This is what you would do in igraph --------------------------------------
#' \dontrun{
#' g   <- ans
#' g[] <- as.character(g)
#' g <- igraph::graph_from_edgelist(g)
#' plot(g, layout = igraph::layout_with_sugiyama(g)[[2]])
#' }
#' 
#' # A performance benchmark with ape::rtree ----------------------------------
#' \dontrun{
#' microbenchmark::microbenchmark(
#' ape = rtree(1e3),
#'   phy = sim_tree(1e3),
#' unit = "relative"
#' )
#' # This is what you would get.
#' Unit: relative
#'   expr     min       lq     mean  median       uq      max neval
#'    ape 14.7598 14.30809 14.30013 16.7217 14.32843 4.754106   100
#'    phy  1.0000  1.00000  1.00000  1.0000  1.00000 1.000000   100
#' }
#' @export
sim_tree <- function(n, edge.length = stats::runif) {
  ans <- .sim_tree(n)
  new_po_tree(
    edges             = ans,
    labels            = attr(ans, "labels"),
    offspring         = attr(ans, "offspring"),
    edge.length       = edge.length(nrow(ans)),
    check.edges       = FALSE,
    check.labels      = FALSE,
    check.offspring   = FALSE,
    check.edge.length = FALSE
    )
}

#' @rdname sim_fun_on_tree
#' @export
sim_fun_on_tree <- function(
  tree,
  psi,
  mu,
  Pi,
  P = 1L
) {
  
  # Must be coerced into a tree of class ape::phylo
  tree <- ape::as.phylo(tree)
  
  # Generating the preorder sequence
  pseq <- ape::postorder(tree)
  pseq <- tree$edge[pseq, 2]
  
  # The preorder is just the inverse of the post order!
  # now, observe that the main function does the call using indexes starting
  # from 0, BUT, that's corrected in the function itself
  pseq <- c(length(tree$tip.label) + 1L, pseq[length(pseq):1L])
  
  # Calling the c++ function that does the hard work
  .sim_fun_on_tree(
    offspring = list_offspring(tree),
    pseq      = pseq,
    psi       = psi,
    mu        = mu,
    Pi        = Pi,
    P         = P
  )
  
}

#' Simulation of Annotated Phylogenetic Trees
#' 
#' @param n Integer scalar. Number of leafs. If not specified, then 
#' @param tree An object of class [=sim_tree::po_tree()].
#' @param P Integer scalar. Number of functions to generate.
#' @template parameters
#' @templateVar psi 1
#' @templateVar mu 1
#' @templateVar Pi 1
#' @return An object of class [=new_aphylo::aphylo()]
#' @family Simulation Functions 
#' @export
#' @examples 
#' # A simple example ----------------------------------------------------------
#' 
#' set.seed(1231)
#' ans <- sim_annotated_tree(n=500)
#' 
sim_annotated_tree <- function(
  n    = NULL,
  tree = NULL,
  P    = 1,
  psi  = c(.05, .05),
  mu   = c(.1,.05),
  Pi   = 1
  ) {
  
  pars <- unlist(c(psi, mu, Pi), recursive = TRUE)
  names(pars) <- c("psi0", "psi1", "mu0", "mu1", "Pi")
  
  # Step 1: Simulate a tree
  if (!length(tree)) {
    
    # Checking if there's n
    if (!length(n))
      stop("When -tree- is not specified, -n- must be specified.")
    
    tree  <- ape::reorder.phylo(ape::rtree(n), "postorder")
    
  } else 
    tree <- as.phylo(tree)
  
  
  # Step 2: Simulate the annotations
  ans <- sim_fun_on_tree(
    tree  = tree,
    psi   = c(pars["psi0"], pars["psi1"]),
    mu    = c(pars["mu0"], pars["mu1"]),
    Pi    = pars["Pi"],
    P     = P
  )
  
  rownames(ans) <- unname(attr(tree, "labels"))
  
  # Creating the aphylo object
  nleaf <- length(tree$tip.label)
  as_aphylo(
    tip.annotation  = ans[1L:nleaf, ,drop=FALSE],
    node.annotation = ans[(nleaf + 1L):nrow(ans), , drop=FALSE],
    tree            = tree
  )
  
}


