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

APHYLO_PARAM_NAMES <- c("psi0", "psi1", "mu0", "mu1", "eta0", "eta1", "Pi")
APHYLO_PARAM_DEFAULT <- structure(
  .Data = c(.1, .1, .1, .1, .9, .9, .1),
  names = APHYLO_PARAM_NAMES
)

#' @rdname aphylo-model
#' @export
aphylo_call <- function() {
  
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
        ans$ll <- ans$ll + log(2^prod(dim(dat$tip.annotation)))
        
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
      names = APHYLO_PARAM_NAMES
      ),
    included = structure(
      .Data = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE),
      names = APHYLO_PARAM_NAMES
    )
  ))
}

#' @rdname aphylo-model
eta <- function(..., env) {
  
  # Adding eta to the objective function
  body(env$fun)[[2]][[3]]$eta  <- bquote(c(p["eta0"], p["eta1"]))
  env$fixed[c("eta0", "eta1")] <- FALSE
  env$included[c("eta0", "eta1")] <- TRUE
  
  # Removing the eta correction
  body(env$fun)[[4]] <- NULL
  
  # Updating
  for (f in unlist(list(...))) {
    body(env$fun)[[2]][[3]]$eta[[2+f]]  <- env$params[paste0("eta", f)]
    env$fixed[paste0("eta", f)] <- TRUE
  }
    
  invisible()
}

validate_dots_in_term <- function(..., expected) {
  
  # Who called me?
  parent_call <- sys.calls()
  parent_call <- parent_call[[length(parent_call) - 1]][[1]]
  parent_call <- as.character(parent_call)
  
  dots <- list(...)
  
  # Capturing dots
  dots <- unlist(list(...))
  if (length(dots) && !all(dots %in% c(0, 1)))
    stop(
      "Argumets passed to `",
      parent_call,
      "`, if any, should be any of {",
      paste0(expected, collapse=", "),
      "}.",
      call. = FALSE
      )
  
  dots
  
}

#' @rdname aphylo-model
psi <- function(..., env) {
  
  # Adding eta to the objective function
  body(env$fun)[[2]][[3]]$psi  <- bquote(c(p["psi0"], p["psi1"]))
  env$fixed[c("psi0", "psi1")] <- FALSE
  env$included[c("psi0", "psi1")] <- TRUE
  
  # Updating
  dots <- validate_dots_in_term(..., expected = c(0,1))
  for (f in dots) {
  
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
  dots <- validate_dots_in_term(..., expected = 0)
  for (f in dots) {
    
    body(env$fun)[[2]][[3]]$Pi <- env$params["Pi"]
    env$fixed[paste0("Pi")] <- TRUE
    
  }
    
  
  invisible()
}

#' @rdname aphylo-model
mu <- function(..., env) {

  # Updating
  dots <- validate_dots_in_term(..., expected = c(0,1))
  for (f in dots) {

    body(env$fun)[[2]][[3]]$mu[[2+f]]  <- env$params[paste0("mu", f)]
    env$fixed[paste0("mu", f)] <- TRUE
    
  }
  
  invisible()
}

#' Basic validations:
#' 1. Are all terms valid terms?
#' 2. Is there any repeated term?
#' @noRd
validate_aphylo_formula <- function(fm) {
  
  # Extracting terms
  fm   <- stats::terms(fm, keep.order = TRUE)
  vars <- attr(fm, "term.labels")
  term_names <- gsub("\\(.+", "", vars)
  
  # Analyzing the terms, are all 
  test <- which(!(term_names %in% unique(gsub("[0-9]", "", APHYLO_PARAM_NAMES))))
  
  if (length(test)) {
    stop("The following terms in the formula are not supported by `aphylo`: '",
         paste0(vars[test], collapse = "', '"), "'.", call. = FALSE)
  }
  
  # Is there any repeated term?
  test <- table(term_names)
  test <- names(test)[which(test > 1)]
  if (length(test)) {
    stop("Some terms are repeated more than once: '",
         paste0(vars[term_names %in% test], collapse = "', '"), "'.",
         call. = FALSE
         )
  }
  

  attr(fm, "variables")
  
  
}

#' This function validates the parameters specified
#' @noRd
validate_parameters <- function(fm, params) {
  
  # Extracting formula terms
  vars <- attr(stats::terms(fm), "term.labels")
  vars <- gsub("\\(.+", "", vars)
  
  # In the case of missing parameters, then return the default according to
  # what the formula is.
  if (missing(params)) {
    
    APHYLO_PARAM_DEFAULT[which(APHYLO_PARAM_NAMES %in% vars)]
    
  } else {
    
    # Checking length
    if (!length(names(params))) {
      
      # If unnamed, are we giving enough information?
      if (length(params) != length(vars))
        stop("The initial parameters have different length than specified in ",
             "the model.", call.=FALSE)
      
      # Matching names by position
      warning("Initial parameres matched by position.")
      
      return(structure(
        .Data = params,
        names = APHYLO_PARAM_NAMES[which(APHYLO_PARAM_NAMES %in% vars)]
      ))
      
    }
    
    # More than needed
    test <- which(!(names(params) %in% vars))
    if (length(test)) {
      stop(
        "The initial parametes must match those specified in the model. ",
        "The following parameters have been overspecified: '",
        paste(names(params)[test], collapse = "', '"), "'. These should ",
        "match the following set of parameters: '",
        paste(vars, collapse = "', '"), "'.", 
        call. = FALSE
        )
    }
    
    # Less than required
    test <- which(!(vars %in% names(params)))
    if (length(test)) {
      stop(
        "The initial parametes must match those specified in the model. ",
        "The following parameters are missing: '",
        paste(vars[test], collapse = "', '"), "'. These should ",
        "match the following set of parameters: '",
        paste(vars, collapse = "', '"), "'.", 
        call. = FALSE
      )
    }
    
    params
    
  }
  
}

#' @rdname aphylo-model
#' @export
aphylo_formula <- function(fm, params) {
  
  # Creating new aphylo call object
  model_call <- aphylo_call()
  
  # Validating formula
  val <- validate_aphylo_formula(fm)
  
  # Validating initial parameters
  params <- validate_parameters(fm, params)
  
  # Is the LHS an aphylo object?
  if (!exists(as.character(val[[2]])))
    stop("The object -",as.character(val[[2]]), "- can't be found.")
    
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
