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
#' aphylo_formula(x ~ mu + psi + Pi)
#' 
#' @name aphylo-model
#' @aliases aphylo-formula
NULL

aphylo_params_names <- c("psi0", "psi1", "mu0", "mu1", "eta0", "eta1", "Pi")

#' @rdname aphylo-model
#' @export
aphylo_call <- function(params) {
  
  if (missing(params))
    params <- structure(rep(0, 7), names = aphylo_params_names)
  
  list2env(
    list(
      fun = function(p, priors, dat, verb_ans = FALSE) {
        
        # Arguments
        args <- list(
          tree = dat,
          psi  = c(0, 0),
          mu   = c(p["mu0"], p["mu1"]),
          eta  = c(.5, .5),
          Pi   = p["mu0"]/(p["mu0"] + p["mu1"])
        )
        
        # Call
        ans <- do.call(
          what = aphylo::LogLike,
          args = args
          )
        
        # Correcting for eta
        ans$ll <- ans$ll*2^prod(dim(dat$tip.annotation))
        
        # Adding priors
        ans$ll <- ans$ll + sum(log(priors(p)))
        
        if (is.infinite(ans$ll))
          ans$ll <- .Machine$double.xmax*sign(ans$ll)*1e-10
        
        # If verbose (not by default)
        if (verb_ans)
          ans
        else
          ans$ll
    },
    fixed    = structure(
      .Data = c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE),
      names = aphylo_params_names
      ),
    included = structure(
      .Data = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE),
      names = aphylo_params_names
    ),
    params   = params
  ))
}

#' @rdname aphylo-model
eta <- function(..., env) {
  
  # Adding eta to the objective function
  body(env$fun)[[2]][[3]]$eta  <- bquote(c(p["eta0"], p["eta1"]))
  env$fixed[c("eta0", "eta1")] <- FALSE
  env$included[c("eta0", "eta1")] <- TRUE
  
  # Removing the eta correction
  body(env$fun)[[4]] <- bquote()
  
  # Updating
  for (f in unlist(list(...))) {
    body(env$fun)[[2]][[3]]$eta[[2+f]]  <- env$params[paste0("eta", f)]
    env$fixed[paste0("eta", f)] <- TRUE
  }
    
  invisible()
}

#' @rdname aphylo-model
psi <- function(..., env) {
  
  # Adding eta to the objective function
  body(env$fun)[[2]][[3]]$psi  <- bquote(c(p["psi0"], p["psi1"]))
  env$fixed[c("psi0", "psi1")] <- FALSE
  env$included[c("psi0", "psi1")] <- TRUE
  
  # Updating
  for (f in unlist(list(...))) {
  
    body(env$fun)[[2]][[3]]$psi[[2+f]]  <- env$params[paste0("psi", f)]
    env$fixed[paste0("psi", f)] <- TRUE
  }
  
  invisible()
}

#' @rdname aphylo-model
Pi <- function(..., env) {
  
  # Adding eta to the objective function
  body(env$fun)[[2]][[3]]$Pi <- bquote(p["Pi"])
  env$fixed["Pi"]            <- FALSE
  
  # Updating (if fixed, then we set whatever value should be included)
  for (f in unlist(list(...))) {
    body(env$fun)[[2]][[3]]$Pi <- env$params["Pi"]
    env$fixed[paste0("Pi")] <- TRUE
  }
    
  
  invisible()
}

#' @rdname aphylo-model
mu <- function(..., env) {
  
  # Updating
  for (f in unlist(list(...))) {
    body(env$fun)[[2]][[3]]$mu[[2+f]]  <- env$params[paste0("mu", f)]
    
    env$fixed[paste0("mu", f)] <- TRUE
  }
  
  invisible()
}

#' @rdname aphylo-model
#' @export
aphylo_formula <- function(fm, params) {
  
  # Creating new aphylo call object
  model_call <- aphylo_call(params)
  
  # Extracting terms
  fm_terms   <- terms(fm, keep.order=TRUE)
  
  # Taking a look at the variables
  val <- attr(fm_terms, "variables")
  
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
  
  # Checking the dimensions
  
  # Returning the model call as a list
  c(
    list(
      model  = fm,
      dat    = eval(val[[2]]),
      params = params
      ),
    as.list(model_call)
  )
    
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
