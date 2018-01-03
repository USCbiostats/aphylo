#' Random tree generation
#' 
#' An alternative to [ape::rtree]. This function was written in C++ and is 
#' significantly faster than `rtree`. 
#' 
#' @param n Integer scalar. Number of leaf nodes.
#' @param edge.length A Function. Used to set the length of the edges.
#' 
#' @details The algorithm was implemented as follows
#' 
#' \enumerate{
#'   \item Initialize `N = {1, ..., n}`, `E` to be empty,
#'   `k = 2*n - 1`
#'   \item While `length(N) != 1` do:
#'   \enumerate{
#'     \item Randomly choose a pair `(i, j)` from `N`
#'     \item Add the edges `E = E U {(k, i), (k, j)}`,
#'     \item Redefine `N = (N \ {i, j}) U {k}`
#'     \item Set `k = k - 1`
#'     \item next
#'   }
#'   \item Use `edge.length(2*n - 1)` (simulating branch lengths).
#' }
#' 
#' 
#' @return An object of class [ape::phylo] with the edgelist as a postorderd,
#' `node.label` and `edge.length`.
#' 
#' @examples
#' # A very simple example ----------------------------------------------------
#' set.seed(1223)
#' newtree <- sim_tree(50)
#' 
#' plot(newtree)
#' 
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
sim_tree <- function(n, edge.length = NULL) {
  
  if (!length(edge.length))
    .sim_tree(n, function(x) rep.int(1, x), FALSE)
  else
    .sim_tree(n, edge.length, TRUE)
  
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
  if (!inherits(tree, "phylo") & !inherits(tree, "aphylo"))
    stop("`tree` must be of class `phylo` or `aphylo`.", call. = FALSE)
  
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
  
  
  # Step 1: Simulate a tree
  if (!length(tree)) {
    
    # Checking if there's n
    if (!length(n))
      stop("When -tree- is not specified, -n- must be specified.")
    tree  <- sim_tree(n)
    
  } else 
    tree <- as.phylo(tree)
  
  
  # Step 2: Simulate the annotations
  ans <- sim_fun_on_tree(
    tree  = tree,
    psi   = psi,
    mu    = mu,
    Pi    = Pi,
    P     = P
  )
  
  # Creating the aphylo object
  nleaf <- length(tree$tip.label)
  as_aphylo(
    tip.annotation  = ans[1L:nleaf, ,drop=FALSE],
    node.annotation = ans[(nleaf + 1L):nrow(ans), , drop=FALSE],
    tree            = tree
  )
  
}


