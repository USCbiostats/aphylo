#' Simulation of Annotated Phylogenetic Trees
#' 
#' @param n Integer scalar. Number of leafs. If not specified, then 
#' @param tree An object of class \code{\link[=sim_tree]{phylo_tree}}.
#' @param P Integer scalar. Number of functions to generate.
#' @param psi A numeric vector of length 2 (see details).
#' @param mu A numeric vector of length 2 (see details).
#' @param Pi A numeric vector of length 2 (see details).
#' (see \code{\link{LogLike}})
#' @return An object of class \code{\link[=get_offspring]{phylo_offspring}}
#' @family Simulation Functions 
#' @export
#' @examples 
#' # A simple example ----------------------------------------------------------
#' 
#' set.seed(1231)
#' ans <- sim_annotated_tree(n=500)
#' 
sim_annotated_tree <- function(
  n = NULL, tree = NULL, P=1,
  psi=c(.05, .05), mu=c(.1,.05), Pi=.5
  ) {
  
  pars <- c(psi0 = psi[1], psi1 = psi[2], mu0 = mu[1], mu1 = mu[2], Pi = Pi)
  
  # Step 1: Simulate a tree
  if (!length(tree)) {
    
    # Checking if there's n
    if (!length(n))
      stop("When -tree- is not specified, -n- must be specified.")
    
    tree <- sim_tree(n)  
    
  } else {
    if (!inherits(tree, "phylo_tree"))
      stop("-tree- must be an object of class -phylo_tree- (see -sim_tree-).")
  }
  
  # Step 2: Simulate the annotations
  ans <- sim_fun_on_tree(
    offspring  = attr(tree, "offspring"),
    noffspring = attr(tree, "noffspring"),
    psi        = c(pars["psi0"], pars["psi1"]),
    mu         = c(pars["mu0"], pars["mu1"]),
    Pi         = pars["Pi"] 
  )
  
  # Creating the phylo_offspring object
  structure(
    list(
      experiment = unname(ans),
      fun_names  = colnames(ans), 
      added      = rep(FALSE, nrow(ans)),
      offspring  = attr(tree, "offspring"),
      noffspring = attr(tree, "noffspring"),
      edgelist   = {
        mat <- array(dim=dim(tree))
        mat[] <- tree[]
        mat
      }
    ),
    class = "phylo_offspring"
  )
}