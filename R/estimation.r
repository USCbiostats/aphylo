
#' Try to compute the inverse of the matrix
#' @return If it fails, then it returns a matrix of size n*n with NAs.
#' @noRd
try_solve <- function(x, ...) {
  
  ans <- tryCatch(solve(x, ...), error = function(e) e)
  
  # If it is an error
  if (inherits(ans, "error")) {
    warning("The algorithm did not converge. Cannot find the inverse of the hessian.")
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
#' @param dat An object of class `new_aphylo` as returned by
#' [new_aphylo()].
#' @param method Character scalar. When `"ABC"`, uses Artificial Bee Colony
#' optimization algorithm, otherwise it uses a method in [stats::optim()]. 
#' @param priors A list of length 3 with functions named `psi`, `mu`,
#' `Pi`
#' @param control A list with parameters for the optimization method (see
#' details).
#' @param lower Numeric vector of length 5. Lower bounds, default to 0.0001.
#' @param upper Numeric vector of length 5. Upper bounds, default to 0.9999
#' @param object An object of class `aphylo_estimates`.
#' @param check.informative Logical scalar. When `TRUE` (default) the algorithm
#' stops with an error when the annotations are uninformative (either 0s or 1s).
#' 
#' @details 
#' 
#' `phylo_mcmc` is a wrapper of [MCMC()], so, instead of treating the
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
#' set.seed(89)
#' dat <- sim_annotated_tree(100, P=2)
#' 
#' # Computing Estimating the parameters 
#' ans  <- aphylo_mle(dat)
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
#' ans_dbeta <- aphylo_mle(dat, priors = mypriors)
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
#'   rep(.1, 7), dat,
#'   control = list(nbatch = 2e5, burnin=1000, thin=200, scale=2e-2)
#' )
#' }
#' 
#' @name aphylo_estimates-class
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
  varcovar
) {
  
  structure(
    list(
      par = par,
      hist = hist,
      ll = ll,
      counts = counts,
      convergence = convergence,
      message = message,
      fun = fun,
      priors = priors,
      dat = dat,
      par0 = par0,
      method = method,
      varcovar = varcovar
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

#' Objective function with priors
#' @noRd
obj_fun_priors <- function(x, dat, priors) {
    
    ll <- LogLike(
      tree       = dat, 
      psi        = x[1:2], 
      mu         = x[3:4], 
      eta        = x[5:6],
      Pi         = x[7],
      verb_ans   = FALSE, 
      check_dims = FALSE
    )$ll + sum(log(priors(x)))
    
    # Checking if we got a infinite result
    if (is.infinite(ll)) return(.Machine$double.xmax*sign(ll)*1e-10)
    ll
}

#' Objective function
#' @noRd
obj_fun <- function(x, dat) {
    
    ll <- LogLike(
      tree       = dat, 
      psi        = x[1:2], 
      mu         = x[3:4],  
      eta        = x[5:6],
      Pi         = x[7],
      verb_ans   = FALSE, 
      check_dims = FALSE
    )$ll
    
    # Checking if we got a infinite result
    if (is.infinite(ll)) return(.Machine$double.xmax*sign(ll)*1e-10)
    ll
}


#' @rdname aphylo_estimates-class
#' @export
aphylo_mle <- function(
  dat,
  method        = "L-BFGS-B",
  priors        = NULL, 
  control       = list(),
  params        = c(rep(.01, 4), .99, .99, .01),
  lower         = 0.0,
  upper         = 1.0,
  check.informative = TRUE
) {
  
  
  # Checking params
  if (!inherits(dat, "aphylo"))
    stop("-dat- should be of class aphylo")
  
  if (length(params) != 7)
    stop("-params- must be of length 5.")
  
  if (!is.numeric(params))
    stop("-params- must be a numeric vector")
  
  # Reducing the peeling sequence
  if (getOption("aphylo_reduce_pseq", FALSE))
    dat$pseq <- reduce_pseq(dat$pseq, with(dat, rbind(tip.annotation, node.annotation)), dat$offspring)

  # If the models is uninformative, then it will return with error
  if (check.informative)
    stop_ifuninformative(dat$tip.annotation)
  
  # In case of fixing parameters
  par0 <- params
  
  # Both available for ABC
  control$fnscale  <- -1
  control$parscale <- rep(1, 7)
  
  # Creating the objective function
  fun <- if (length(priors)) {
    function(x, dat) {

      ll <- LogLike(
        tree       = dat, 
        psi        = x[1:2], 
        mu         = x[3:4], 
        eta        = x[5:6],
        Pi         = x[7],
        verb_ans   = FALSE, 
        check_dims = FALSE
        )$ll + sum(log(priors(x)))
      
      # Checking if we got a infinite result
      if (is.infinite(ll)) return(.Machine$double.xmax*sign(ll)*1e-10)
      ll
    }
  } else {
    function(x, dat) {

      ll <- LogLike(
        tree       = dat, 
        psi        = x[1:2], 
        mu         = x[3:4],  
        eta        = x[5:6],
        Pi         = x[7],
        verb_ans   = FALSE, 
        check_dims = FALSE
      )$ll

      # Checking if we got a infinite result
      if (is.infinite(ll)) return(.Machine$double.xmax*sign(ll)*1e-10)
      ll
    }
  }
  
  # Optimizing
  ans <- do.call(
    stats::optim, 
    c(
      list(par = params, fn = fun, method=method, upper = upper, lower = lower,
           hessian=FALSE, dat = dat, control=control)
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
  hessian <- stats::optimHess(ans$par, fun, dat = dat, control = control)
  
  # Working on answer
  env <- new.env()
  environment(fun)    <- env
  environment(dat)    <- env
  environment(par0)   <- env
  if (length(priors)) 
    environment(priors) <- env
  
  # Naming the parameter estimates
  names(ans$par) <- c("psi0", "psi1", "mu0", "mu1", "eta0", "eta1", "Pi")
  
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
    fun        = fun,
    priors     = priors,
    dat        = dat,
    par0       = par0,
    method     = method,
    varcovar   = try_solve(-hessian, tol = 1e-100)
  )
}

#' @export
#' @rdname aphylo_estimates-class
print.aphylo_estimates <- function(x, ...) {
  # Function to print a bar with variable width
  catbar <- function() paste0(rep("-",options()$width), collapse="")
  
  sderrors   <- sqrt(diag(x$varcovar))

  with(x, {
    cat(
      sep = "\n",
      "\nESTIMATION OF ANNOTATED PHYLOGENETIC TREE",
      sprintf(
        "ll: %9.4f,\nMethod used: %s (%i iterations)", ll, method, x$counts),
      if (method == "mcmc") 
        NULL
      else
        sprintf("convergence: %i (see ?optim)", convergence)
      ,
      sprintf("Leafs\n # of Functions %i", ncol(dat$tip.annotation)),
      #paste0(sprintf(" # of %s: %5i (%2.0f%%)", names(props), props, propspcent), collapse="\n"),
            "\n         Estimate  Std. Error",
      sprintf(" psi[0]    %6.4f      %6.4f", par["psi0"], sderrors["psi0"]),
      sprintf(" psi[1]    %6.4f      %6.4f", par["psi1"], sderrors["psi1"]),
      sprintf(" mu[0]     %6.4f      %6.4f", par["mu0"], sderrors["mu0"]),
      sprintf(" mu[1]     %6.4f      %6.4f", par["mu1"], sderrors["mu1"]),
      sprintf(" eta[0]    %6.4f      %6.4f", par["eta0"], sderrors["eta0"]),
      sprintf(" eta[1]    %6.4f      %6.4f", par["eta1"], sderrors["eta1"]),
      sprintf(" Pi        %6.4f      %6.4f\n", par["Pi"], sderrors["Pi"])
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

#' @param x An object of class `aphylo_estimates`.
#' @param ... Further arguments passed to the method.
#' @rdname aphylo_estimates-class
#' @export
plot.aphylo_estimates <- function(
  x,
  ...
  ) {
  
  plot_LogLike.aphylo_estimates(x, ...)
}

#' @rdname aphylo_estimates-class
#' @return In the case of `aphylo_mcmc`, `hist` is an object of class
#' [coda::mcmc.list()].
#' @export
aphylo_mcmc <- function(
  params,
  dat,
  priors        = NULL,
  control       = list(),
  check.informative = TRUE
) {
  
  # Checking control
  if (!length(control$nbatch)) control$nbatch <- 2e3
  if (!length(control$scale))  control$scale  <- .01
  if (!length(control$ub))     control$ub     <- rep(1, 7)
  if (!length(control$lb))     control$lb     <- rep(0, 7)
  if (!length(control$useCpp)) control$useCpp <- TRUE
  
  # Checking params
  if (!inherits(dat, "aphylo"))
    stop("-dat- should be of class aphylo")
  
  if (length(params) != 7)
    stop("-params- must be of length 5.")
  
  if (!is.numeric(params))
    stop("-params- must be a numeric vector")
  
  # If the models is uninformative, then it will return with error
  if (check.informative)
    stop_ifuninformative(dat$tip.annotation)
  
  # In case of fixing parameters
  par0 <- params
  
  # Creating the objective function
  fun <- if (length(priors)) {
    function(params, dat) {

      LogLike(
        tree       = dat,
        psi        = params[1:2] ,
        mu         = params[3:4] ,
        eta        = params[5:6] ,
        Pi         = params[7],
        verb_ans   = FALSE, 
        check_dims = FALSE
      )$ll + sum(log(priors(params)))
      
    }
  } else {
    function(params, dat) {

      LogLike(
        tree       = dat,
        psi        = params[1:2] ,
        mu         = params[3:4] ,
        eta        = params[5:6] ,
        Pi         = params[7],
        verb_ans   = FALSE, 
        check_dims = FALSE
      )$ll
      
    }
  }
  
  # Naming the parameter estimates
  names(params) <- c("psi0", "psi1", "mu0", "mu1", "eta0", "eta1", "Pi")
  
  # Running the MCMC
  ans <- do.call(amcmc::MCMC, c(list(fun = fun, initial = params, dat=dat), control))
  
  # We treat all chains as mcmc.list
  if (!inherits(ans, "mcmc.list"))
    ans <- coda::mcmc.list(ans)
  
  # Working on answer
  env <- new.env()
  environment(fun)    <- env
  environment(dat)    <- env
  environment(par0)   <- env
  if (length(priors)) 
    environment(priors) <- env
  
  # Returning
  new_aphylo_estimates(
    par        = colMeans(do.call(rbind, ans)),
    hist       = ans,
    ll         = mean(apply(do.call(rbind, ans), 1, fun, dat=dat), na.rm=TRUE),
    counts     = control$nbatch,
    convergence = NA,
    message    = NA,
    fun        = fun,
    priors     = priors,
    dat        = dat,
    par0       = par0,
    method     = "mcmc",
    varcovar   = var(do.call(rbind, ans))
  )
}
