#' Formulas in `aphylo`
#' 
#' @param ... Either 0, 1 or both. Depending on the parameter, the index of the
#' model parameter that will be set as fixed.
#' @param env Environment (not to be called by the user).
#' @return A list with the following elements:
#' * `fun` A function. The log-likelihood function.
#' * `fixed` Logical vector. 
#' 
#' @examples
#' set.seed(12)
#' x <- sim_annotated_tree(10)
#' 
#' # Baseline model
#' aphylo_formula(x ~ mu)
#' 
#' # Mislabeling probabilities
#' aphylo_formula(x ~ mu + psi)
#' 
#' # Mislabeling probabilities and etas(fixed)
#' aphylo_formula(x ~ mu + psi + eta(0, 1))
#' 
#' # Mislabeling probabilities and Pi 
#' #aphylo_formula(x ~ mu + psi + Pi)
#' 
#' @name aphylo-model
#' @aliases aphylo-formula
NULL

#' @rdname aphylo-model
#' @export
aphylo_call <- function() {
  list2env(
    list(
      fun = function(p, dat) {
      
      # Arguments
      args <- list(
        tree = dat,
        psi  = c(0, 0),
        mu   = p[c("mu0", "mu1")],
        eta  = c(.5, .5),
        Pi   = p["mu0"]/(p["mu0"] + p["mu1"])
      )
      
      # Call
      do.call(
        what = aphylo::LogLike,
        args = args
        )$ll + priors(p)
    },
    fixed = structure(
      .Data = c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE),
      names = c("psi0", "psi1", "mu0", "mu1", "eta0", "eta1", "Pi")
      )
  ))
}

#' @rdname aphylo-model
eta <- function(..., env) {
  
  # Adding eta to the objective function
  body(env$fun)[[2]][[3]]$eta  <- bquote(p[c("eta0", "eta1")])
  env$fixed[c("eta0", "eta1")] <- FALSE
  
  # Updating
  for (f in unlist(list(...))) 
    env$fixed[paste0("eta", f)] <- TRUE
    
  invisible()
}

#' @rdname aphylo-model
psi <- function(..., env) {
  
  # Adding eta to the objective function
  body(env$fun)[[2]][[3]]$psi  <- bquote(p[c("psi0", "psi1")])
  env$fixed[c("psi0", "psi1")] <- FALSE
  
  # Updating
  for (f in unlist(list(...))) 
    env$fixed[paste0("psi", f)] <- TRUE
  
  invisible()
}

#' @rdname aphylo-model
Pi <- function(..., env) {
  
  # Adding eta to the objective function
  body(env$fun)[[2]][[3]]$Pi <- bquote(p["Pi"])
  env$fixed["Pi"]            <- FALSE
  
  # Updating
  for (f in unlist(list(...))) 
    env$fixed[paste0("Pi")] <- TRUE
  
  invisible()
}

#' @rdname aphylo-model
mu <- function(..., env) {
  
  # Updating
  for (f in unlist(list(...))) 
    env$fixed[paste0("mu", f)] <- TRUE
  
  invisible()
}

#' @rdname aphylo-model
#' @export
aphylo_formula <- function(fm) {
  
  # Creating new aphylo call object
  model_call <- aphylo_call()
  
  # Extracting terms
  fm <- terms(fm, keep.order=TRUE)
  
  # Taking a look at the variables
  val <- attr(fm, "variables")
  
  # Is the LHS an aphylo object?
  if (!inherits(eval(val[[2]]), "aphylo"))
    stop("The LHS of the equation should be an `aphylo` object.", call. = FALSE)
  
  # Mofiying the likelihood function and the parameters for the mcmc
  for (i in 3:length(val))
    if (!is.call(val[[i]])) {
      eval(
        call(as.character(val[[i]]), env = bquote(model_call))
      )
    } else {
      val[[i]]$env <- bquote(model_call)
      eval(val[[i]])
    }
  
  # Returning the model call as a list
  as.list(model_call)
}



# 
# e <- list2env(list(x = x))
# x <- 1
# 
# z <- function() {
#   summary(x)
# }
# 
# environment(z) <- e
# z()
