#' Formulas in `aphylo`
#' 
#' @param ... Either 0, 1 or both. Depending on the parameter, the index of the
#' model parameter that will be set as fixed.
#' @param fm A formula. Model of the type `<aphylo-object> ~ <parameters>` (see 
#' examples).
#' @param priors (optional) A function. Prior for the model.
#' @param env Environment (not to be called by the user).
#' @param params Numeric vector with model parameters.
#' @return A list with the following elements:
#' * `fun` A function. The log-likelihood function.
#' * `fixed` Logical vector. 
#' 
#' @examples
#' set.seed(12)
#' x <- raphylo(10)
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
aphylo_call <- function(params, priors) {
  
  ans <- list2env(
    list(
      fun = function(p, dat, priors, verb_ans = FALSE) {
        
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
        ans$ll <- ans$ll +
          0.69314718055994528623*sum(Nann(dat))*sum(Nannotated(dat))
        
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
    fixed = structure(
      .Data = rep(FALSE, length(params)),
      names = names(params)
      ),
    params = params
  ))
  
  if (!missing(priors))
    formals(ans$fun)$priors <- priors
  
  ans
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
      "Arguments passed to `",
      parent_call,
      "`, if any, should be any of {",
      paste0(expected, collapse=", "),
      "}.",
      call. = FALSE
      )
  
  dots
  
}

#' @rdname aphylo-model
eta <- function(..., env) {
  
  # Adding eta to the objective function
  body(env$fun)[[2]][[3]]$eta  <- bquote(c(p["eta0"], p["eta1"]))
  env$fixed[c("eta0", "eta1")] <- FALSE
  
  # Removing the eta correction
  body(env$fun)[[4]] <- NULL
  
  # Updating
  dots <- validate_dots_in_term(..., expected = c(0,1))
  for (f in dots) 
    env$fixed[paste0("eta", f)] <- TRUE
  
  invisible()
}

#' @rdname aphylo-model
psi <- function(..., env) {
  
  # Adding eta to the objective function
  body(env$fun)[[2]][[3]]$psi  <- bquote(c(p["psi0"], p["psi1"]))
  env$fixed[c("psi0", "psi1")] <- FALSE

  # Updating
  dots <- validate_dots_in_term(..., expected = c(0,1))
  for (f in dots) 
    env$fixed[paste0("psi", f)] <- TRUE
  
  invisible()
}

#' @rdname aphylo-model
Pi <- function(..., env) {
  
  # Adding eta to the objective function
  body(env$fun)[[2]][[3]]$Pi <- bquote(p["Pi"])
  env$fixed["Pi"]            <- FALSE
  
  # Updating (if fixed, then we set whatever value should be included)
  dots <- validate_dots_in_term(..., expected = 0)
  for (f in dots)
    env$fixed[paste0("Pi")] <- TRUE
  
  invisible()
}

#' @rdname aphylo-model
mu <- function(..., env) {

  # Updating
  dots <- validate_dots_in_term(..., expected = c(0,1))
  for (f in dots)
    env$fixed[paste0("mu", f)] <- TRUE
  
  invisible()
}

#' Basic validations:
#' 1. Are all terms valid terms?
#' 2. Is there any repeated term?
#' @noRd
validate_aphylo_formula <- function(fm) {
  
  # Extracting terms
  fm1  <- stats::terms(fm, keep.order = TRUE)
  vars <- attr(fm1, "term.labels")
  term_names <- gsub("\\(.+", "", vars)
  
  # Analyzing the terms, are all 
  TNAMES <- unique(gsub("[0-9]", "", APHYLO_PARAM_NAMES))
  test. <- which(!(term_names %in% TNAMES))
  
  if (length(test.)) {
    stop("The following terms in the formula are not supported by `aphylo`: '",
         paste0(vars[test.], collapse = "', '"), "'. ", 
         "Only `", paste(TNAMES, collapse="()`, `"), "()` are allowed.",
         call. = FALSE)
  }
  
  # Is there any repeated term?
  test. <- table(term_names)
  test. <- names(test.)[which(test. > 1)]
  if (length(test.)) {
    stop("Some terms are repeated more than once: '",
         paste0(vars[term_names %in% test.], collapse = "', '"), "'.",
         call. = FALSE
         )
  }
  
  # Is mu present?
  if (!("mu" %in% term_names))
    fm <- stats::update.formula(fm, ~. + mu)
  
  fm
  
  
}

#' This function validates the parameters specified
#' @noRd
validate_parameters <- function(fm, params) {
  
  # Extracting formula terms
  vars <- attr(stats::terms(fm), "term.labels")
  vars <- gsub("\\(.+", "", vars)
  
  
  # Getting the default
  vars <- paste0("^(",paste0(vars, collapse="|"), ")")
  vars <- APHYLO_PARAM_NAMES[grepl(vars, x = APHYLO_PARAM_NAMES)]
  
  # In the case of missing parameters, then return the default according to
  # what the formula is.
  if (missing(params)) {
    
    message("No parameters were specified. Default will be used instead.")
    APHYLO_PARAM_DEFAULT[vars]
    
  } else {
    
    # Retrieving name of the parameters
    pnames <- switch(
      class(params),
      matrix = colnames(params),
      names(params)
      )
    
    # Checking length
    if (!length(pnames)) {
      
      # If unnamed, are we giving enough information?
      npars <- switch (class(params),
        matrix = ncol(params),
        length(params)
      )
      if (npars != length(vars))
        stop("The initial parameters have different length than specified in ",
             "the model.", call.=FALSE)
      
      # Matching names by position
      warning("Initial parameres matched by position.", call. = FALSE)
      
      if (is.matrix(params)) 
        return(structure(.Data = params, dimnames = list(NULL, vars)))
      else
        return(structure(.Data = params, names = vars))
      
    }
    
    # More than needed
    test. <- which(!(pnames %in% vars))
    if (length(test.)) {
      stop(
        "The initial parametes must match those specified in the model. ",
        "The following parameters have been overspecified: '",
        paste(pnames[test.], collapse = "', '"), "'. These should ",
        "match the following set of parameters: '",
        paste(vars, collapse = "', '"), "'.", 
        call. = FALSE
        )
    }
    
    # Less than required
    test. <- which(!(vars %in% pnames))
    if (length(test.)) {
      stop(
        "The initial parametes must match those specified in the model. ",
        "The following parameters are missing: '",
        paste(vars[test.], collapse = "', '"), "'. These should ",
        "match the following set of parameters: '",
        paste(vars, collapse = "', '"), "'.", 
        call. = FALSE
      )
    }
    
    # Sorting accordingly
    if (is.matrix(params))
      params[,vars,drop=FALSE]
    else
      params[vars]
    
  }
  
}

#' @rdname aphylo-model
#' @export
aphylo_formula <- function(fm, params, priors, env = parent.frame()) {
  
  # Validating formula
  fm  <- validate_aphylo_formula(fm)
  val <- attr(stats::terms(fm), "variables")
  
  # Validating initial parameters
  params <- validate_parameters(fm, params)
  
  # Creating new aphylo call object
  model_call <- aphylo_call(params, priors)
  
  # Is the LHS an aphylo object?
  if (!exists(as.character(val[[2]]), envir = env))
    stop("The object -", as.character(val[[2]]), "- can't be found.", call. = FALSE)
  
  LHS <- eval(val[[2]], envir = env)
  if (!is.aphylo(LHS) && !is.multiAphylo(LHS))
    stop(
      "The LHS of the equation should be either a list or a single ",
      "aphylo object.", call. = FALSE)
  
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
  c(
    list(
      model  = fm,
      dat    = eval(val[[2]], env)
      ),
    as.list(model_call)
  )
    
}
