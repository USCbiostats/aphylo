
#' Try to compute the inverse of the matrix
#' @return If it fails, then it returns a matrix of size n*n with NAs.
#' @noRd
try_solve <- function(x, ...) {
  
  ans <- tryCatch(MASS::ginv(x, ...), error = function(e) e)
  
  # If it is an error
  if (inherits(ans, "error")) {
    warning("The algorithm did not converge. Cannot find the inverse of the hessian.",
            call. = FALSE)
    return(matrix(ncol = length(x), nrow=length(x)))
  }
    
  
  ans
}

#' Parameter estimation of Annotated Phylogenetic Trees
#'
#' The optimization is done via `optim`.
#'
#' @param params A vector of length 7 with initial parameters. In particular
#' `psi[1]`, `psi[2]`, `mu[1]`, `mu[2]`, `eta[1]`, `eta[2]` and `Pi`.
#' @param model A model as specified in [aphylo-model].
#' @param method Character scalar. When `"ABC"`, uses Artificial Bee Colony
#' optimization algorithm, otherwise it uses a method in [stats::optim()]. 
#' @param priors A function to be used as prior for the model (see [bprior]).
#' @param control A list with parameters for the optimization method (see
#' details).
#' @param lower,upper Numeric vectors defining the lower and upper bounds respectively.
#' @param object,x An object of class `aphylo_estimates`.
#' @param check.informative Logical scalar. When `TRUE` the algorithm
#' stops with an error when the annotations are uninformative (either 0s or 1s).
#' @param reduced_pseq Logical. When `TRUE` it will use a reduced peeling sequence
#' in which it drops unannotated leafs. If the model includes `eta` this is set
#' to `FALSE`.
#' @param ... Further arguments passed to the method
#' 
#' @details 
#' 
#' `phylo_mcmc` is a wrapper of [amcmc::MCMC], so, instead of treating the
#' problem as a maximization problem, `phylo_mcmc` generates a **Markov Chain**.
#' The default values of `control` are:
#' 
#' \tabular{ll}{
#' `nbatch` \tab Integer scalar. Number of mcmc steps. Default `2e3`. \cr
#' `scale` \tab Numeric scalar. Default `0.01`. \cr
#' `lb` \tab Numeric vector. Default `rep(1e-20, 5)`. \cr
#' `ub` \tab Numeric vector. Default `rep(1 - 1e-20, 5)`. \cr
#' }
#' 
#' @return 
#' A list of class `aphylo_estimates` with the following elements:
#' \item{par}{A numeric vector of length 5 with the solution.}
#' \item{hist}{A numeric matrix of size `counts*5` with the solution path (length 2 if used `optim`
#' as the intermediate steps are not available to the user).}
#' \item{ll}{A numeric scalar with the value of `fun(par, dat)`. The value of the log likelihood.}
#' \item{counts}{Integer scalar number of steps/batch performed.}
#' \item{convergence}{Integer scalar. Equal to 0 if `optim` converged. See `optim`.}
#' \item{message}{Character scalar. See `optim`.}
#' \item{fun}{A function (the objective function).}
#' \item{priors}{If specified, the function `priors` passed to the method.}
#' \item{dat}{The data `dat` provided to the function.}
#' \item{par0}{A numeric vector of length 5 with the initial parameters.}
#' \item{method}{Character scalar with the name of the method used.}
#' \item{varcovar}{A matrix of size 5*5. The estimated covariance matrix.}
#' 
#' @examples 
#' 
#' # Using simulated data ------------------------------------------------------
#' set.seed(19)
#' dat <- sim_annotated_tree(100)
#' dat <- rdrop_annotations(dat, .4)
#' 
#' # Computing Estimating the parameters 
#' ans  <- aphylo_mle(dat ~ psi + mu + eta + Pi)
#' ans
#' 
#' # Plotting the path
#' plot(ans)
#' 
#' # Computing Estimating the parameters Using Priors for all the parameters
#' mypriors <- function(params) {
#'     dbeta(params, c(2, 2, 2, 2, 1, 10, 2), rep(10, 7))
#' }
#' 
#' ans_dbeta <- aphylo_mle(dat ~ psi + mu + eta + Pi, priors = mypriors)
#' ans_dbeta
#' 
#' 
#' # Using the MCMC ------------------------------------------------------------
#' 
#' \dontrun{
#' 
#' set.seed(1233)
#' # Simulating a tree
#' tree <- sim_tree(200)
#' 
#' # Simulating functions
#' dat <- sim_annotated_tree(
#'   tree = tree,
#'   psi  = c(.01, .03),
#'   mu   = c(.05, .02),
#'   Pi   = .5
#' )
#' 
#' # Running the MCMC
#' set.seed(1231)
#' 
#' ans_mcmc <- aphylo_mcmc(
#'   dat ~ mu + psi + eta + Pi,
#'   control = list(nbatch = 2e5, burnin=1000, thin=200, scale=2e-2)
#' )
#' }
#' 
#' @name aphylo_estimates-class
#' @aliases aphylo_estimates
NULL

new_aphylo_estimates <- function(
  par,
  hist,
  ll,
  counts,
  convergence,
  message,
  fun,
  priors,
  dat,
  par0,
  method,
  varcovar,
  call
) {
  
  structure(
    list(
      par         = par,
      hist        = hist,
      ll          = ll,
      counts      = counts,
      convergence = convergence,
      message     = message,
      fun         = fun,
      priors      = priors,
      dat         = dat,
      par0        = par0,
      method      = method,
      varcovar    = varcovar,
      call        = call
    ),
    class = "aphylo_estimates"
  )
}

#' Stop if the model is uninformative
#' @noRd
stop_ifuninformative <- function(tip.annotation) {
  tab <- fast_table_using_labels(tip.annotation, c(0L, 1L))
  if (tab[1L] == 0L | tab[2] == 0L)
    stop("The model is uninformative (either there's only 0s or 1s).", call. = FALSE)
}

#' @rdname aphylo_estimates-class
#' @export
aphylo_mle <- function(
  model,
  params,
  method            = "L-BFGS-B",
  priors            = function(p) 1, 
  control           = list(),
  lower             = 1e-10,
  upper             = 1 - 1e-10,
  check.informative = getOption("aphylo.informative", FALSE),
  reduced_pseq      = getOption("aphylo_reduce_pseq", FALSE)
) {
  
  # Getting the call
  cl <- match.call()
  
  # Parsing the formula
  env   <- parent.frame()
  model <- aphylo_formula(model, params, priors, env = env)
  
  if (any(model$fixed))
    warnings("Fixing parameters is ignored in MLE estimation.")
  
  # Reducing the peeling sequence
  # This only happens if the eta parameter is not included
  if (length(model$fixed["eta0"]))
    reduced_pseq <- FALSE
  
  if (reduced_pseq)
    model$dat$pseq <- reduce_pseq(
      model$dat$pseq,
      with(model$dat, rbind(tip.annotation, node.annotation)),
      model$dat$offspring
    )
  
  # If the models is uninformative, then it will return with error
  if (check.informative)
    stop_ifuninformative(model$dat$tip.annotation)
  
  # Both available for ABC
  control$fnscale  <- -1
  control$parscale <- rep(1, length(model$params))
  
  # Optimizing
  ans <- do.call(
    stats::optim, 
    c(
      list(
        par      = model$params,
        fn       = model$fun,
        dat      = model$dat,
        priors   = priors,
        verb_ans = FALSE,
        method   = method,
        upper    = upper,
        lower    = lower,
        hessian  = FALSE,
        control = control
        )
      )
    )
  
  ans <- list(
    par         = ans$par,
    value       = ans$value,
    convergence = ans$convergence,
    message     = ans$message,
    counts      = ans$counts["function"]
    )
    
    
  # Computing the hessian (information matrix)
  hessian <- stats::optimHess(
    ans$par, model$fun, dat = model$dat, priors = priors, verb_ans = FALSE,
    control = control
    )
  
  # Hessian for observed information matrix
  dimnames(hessian) <- list(names(ans$par), names(ans$par))
  
  # Returning
  new_aphylo_estimates(
    par        = ans$par,
    hist       = NULL,
    ll         = ans$value,
    counts     = ans$counts,
    convergence = ans$convergence,
    message    = ans$message,
    fun        = model$fun,
    priors     = priors,
    dat        = model$dat,
    par0       = model$params,
    method     = method,
    varcovar   = try_solve(-hessian, tol = 1e-100),
    call       = cl
  )
}

#' @export
#' @rdname aphylo_estimates-class
print.aphylo_estimates <- function(x, ...) {
  # Function to print a bar with variable width
  catbar <- function() paste0(rep("-",options()$width), collapse="")
  
  sderrors   <- sqrt(diag(x$varcovar))
  
  ans <- sprintf("\n Leafs:\n # of Functions %i", ncol(x$dat$tip.annotation))
  ans <- c(ans, sprintf("\n %-6s  %6s  %6s", "", "Estimate", "Std. Err."))
  for (p in names(x$par)) {
    ans <- c(
      ans,
      with(x, sprintf("\n %-6s  %6.4f    %6.4f", p, par[p], sderrors[p]))
      )
  }

  with(x, {
    cat(
      sep = "",
      "\nESTIMATION OF ANNOTATED PHYLOGENETIC TREE\n",
      "\n Call: ", paste(deparse(x$call), sep="\n", collapse="\n"), 
      sprintf(
        "\n ll: %-9.4f,\n Method used: %s (%i iterations)", ll, method, x$counts),
      if (method == "mcmc") 
        NULL
      else
        sprintf("\n convergence: %i (see ?optim)", convergence)
      ,
      ans, "\n\n"
      )
    
  })
  
  invisible(x)
}

#' @export
#' @rdname aphylo_estimates-class
coef.aphylo_estimates <- function(object, ...) {
  object$par
}

#' @export
#' @rdname aphylo_estimates-class
vcov.aphylo_estimates <- function(object, ...) {
  object$varcovar
}

#' @rdname plot_logLik
#' @export
#' @include plot_logLik.r
plot.aphylo_estimates <- plot_logLik.aphylo_estimates

#' @export
logLik.aphylo_estimates <- function(object, ...) {
  
  ans <- with(object, fun(par, priors = priors, dat = dat, verb_ans = TRUE))
  
  structure(
    .Data = ans$ll,
    class = "logLik",
    df    = length(object$par),
    Pr    = ans$Pr
  )
  
}

APHYLO_DEFAULT_MCMC_CONTROL <- list(
  nbatch    = 1e5L,
  burnin    = 1e4L,
  thin      = 20L,
  scale     = .01,
  ub        = 1,
  lb        = 0,
  nchains   = 2L,
  multicore = TRUE,
  useCpp    = FALSE,
  autostop  = 5e3
)

#' @rdname aphylo_estimates-class
#' @return In the case of `aphylo_mcmc`, `hist` is an object of class
#' [coda::mcmc.list()].
#' @export
aphylo_mcmc <- function(
  model,
  params,
  priors            = uprior(),
  control           = list(),
  check.informative = getOption("aphylo.informative", FALSE),
  reduced_pseq      = getOption("aphylo_reduce_pseq", FALSE)
) {
  
  # Getting the call
  cl <- match.call()
  
  # Parsing the formula
  env   <- parent.frame()
  model <- aphylo_formula(model, params, priors, env = env)

  # Reducing the peeling sequence
  # This only happens if the eta parameter is not included
  if ("eta0" %in% names(model$fixed))
    reduced_pseq <- FALSE
    
  if (reduced_pseq)
    model$dat$pseq <- reduce_pseq(
      model$dat$pseq,
      with(model$dat, rbind(tip.annotation, node.annotation)),
      model$dat$offspring
    )
  
  # Checking control
  for (n in names(APHYLO_DEFAULT_MCMC_CONTROL)) {
    if (!length(control[[n]]))
      control[[n]] <- APHYLO_DEFAULT_MCMC_CONTROL[[n]]
  }
  
  # Fixed parameters
  if (!length(control$fixed))  control$fixed  <- model$fixed

  # If the models is uninformative, then it will return with error
  if (check.informative)
    stop_ifuninformative(model$dat$tip.annotation)
  
  # Running the MCMC
  ans <- do.call(
    amcmc::MCMC, 
    c(
      list(
        fun      = model$fun,
        dat      = model$dat,
        priors   = priors,
        verb_ans = FALSE,
        initial  = model$params
        ),
      control
      )
    )
  
  # We treat all chains as mcmc.list
  if (!inherits(ans, "mcmc.list"))
    ans <- coda::mcmc.list(ans)
  
  # Working on answer
  par <- colMeans(do.call(rbind, ans))
  
  # Returning
  new_aphylo_estimates(
    par         = par,
    hist        = ans,
    ll          = do.call(model$fun, list(dat = model$dat, priors=priors, p=par, verb_ans=FALSE)),
    counts      = coda::niter(ans),
    convergence = NA,
    message     = NA,
    fun         = model$fun,
    priors      = priors,
    dat         = model$dat,
    par0        = model$params,
    method      = "mcmc",
    varcovar    = var(do.call(rbind, ans)),
    call        = cl
  )
}
