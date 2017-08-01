#' Maximum Likelihood Estimation of the Phylogenetic Tree
#'
#' The optimization is done either via \code{abc_cpp} or \code{optim}.
#'
#' @param params A vector of length 5 with initial parameters. In particular
#' \code{psi[1]}, \code{psi[2]}, \code{mu[1]}, \code{mu[2]}, and \code{Pi}.
#' @param dat An object of class \code{new_aphylo} as returned by
#' \code{\link{new_aphylo}}.
#' @param method Character scalar. When \code{"ABC"}, uses Artificial Bee Colony
#' optimization algorithm, otherwise it uses a method in \code{\link[stats:optim]{optim}}. 
#' @param priors A list of length 3 with functions named \code{psi}, \code{mu},
#' \code{Pi}
#' @param fix.params A Logical vector of length 5. Whether or not to fix
#' a particular parameter to that of what was specified in \code{params}.
#' @param control A list with parameters for the optimization method (see
#' details).
#' @param lower Numeric vector of length 5. Lower bounds, default to 0.0001.
#' @param upper Numeric vector of length 5. Upper bounds, default to 0.9999
#' @param object An object of class \code{phylo_mle}.
#' 
#' @details When \code{method="ABC"}, the optimization is done via 
#' \bold{Artificial Bee Colony method}. The default \code{control} parameters, which
#' are passed to \code{\link[ABCoptim:abc_cpp]{abc_cpp}} in the \pkg{ABCoptim} package,
#' are the following:
#' 
#' \tabular{ll}{
#' \code{criter} \tab Integer scalar. Default \code{50L}. \cr
#' \code{maxCycle} \tab Integer scalar. Default \code{500L}. \cr
#' \code{lb} \tab Numeric scalar. Default \code{1e-20}. \cr
#' \code{ub} \tab Numeric scalar. Default \code{1 - 1e-20}.
#' }
#' 
#' In the case of \code{method != "ABC"}, the algorithm is somewhat numerically
#' unstable when using starting points with values higher than \code{0.1}.
#' It is recomended to use values closer to \code{0.1} instead (as the default).
#' 
#' \code{phylo_mcmc} is a wrapper of \code{\link{MCMC}}, so, instead of treating the
#' problem as a maximization problem, \code{phylo_mcmc} generates a \bold{Markov Chain}.
#' The default values of \code{control} are:
#' 
#' \tabular{ll}{
#' \code{nbatch} \tab Integer scalar. Number of mcmc steps. Default \code{2e3}. \cr
#' \code{scale} \tab Numeric scalar. Default \code{0.01}. \cr
#' \code{lb} \tab Numeric vector. Default \code{rep(1e-20, 5)}. \cr
#' \code{ub} \tab Numeric vector. Default \code{rep(1 - 1e-20, 5)}. \cr
#' }
#' 
#' @return 
#' A list of class \code{phylo_mle} with the following elements:
#' \item{par}{A numeric vector of length 5 with the solution.}
#' \item{hist}{A numeric matrix of size \code{counts*5} with the solution path (length 2 if used \code{optim}
#' as the intermediate steps are not available to the user).}
#' \item{ll}{A numeric scalar with the value of \code{fun(par, dat)}. The value of the log likelihood.}
#' \item{counts}{Integer scalar number of steps/batch performed.}
#' \item{convergence}{Integer scalar. Equal to 0 if \code{optim} converged. See \code{optim}.}
#' \item{message}{Character scalar. See \code{optim}.}
#' \item{fun}{A function (the objective function).}
#' \item{priors}{If specified, the function \code{priors} passed to the method.}
#' \item{dat}{The data \code{dat} provided to the function.}
#' \item{par0}{A numeric vector of length 5 with the initial parameters.}
#' \item{fix.params}{Logical vector of length 5, as passed to the method.}
#' \item{method}{Character scalar with the name of the method used.}
#' \item{varcovar}{A matrix of size 5*5. The estimated covariance matrix.}
#' 
#' @examples 
#' 
#' # Using simulated data ------------------------------------------------------
#' set.seed(890)
#' dat <- sim_annotated_tree(100, P=2)
#' 
#' 
#' # Computing Estimating the parameters 
#' # for some reason, starting with parameters equal to .5 breaks NR.
#' ans_nr  <- phylo_mle(dat)
#' ans_abc <- phylo_mle(dat, method = "ABC")
#' 
#' # Plotting the path
#' with(ans_nr, plot(
#'  - apply(hist, 1, fun, dat=dat),
#'  type = "l",
#'  xlab = "Step",
#'  ylab = "Log-Likelihood"
#' ))
#' 
#' 
#' # Computing Estimating the parameters Using Priors for PSI 
#' mypriors <- function(params) {
#'     dbeta(params[1:2], 2, 10)
#' }
#' ans_nr_dbeta <- phylo_mle(dat, priors = mypriors)
#' ans_abc_dbeta <- phylo_mle(dat, priors = mypriors, method = "ABC")
#' 
#' # Plotting the path
#' oldpar <- par(no.readonly = TRUE)
#' par(mfrow = c(2, 2))
#' 
#' plot(ans_nr, main = "No priors L-BFGS-B")
#' plot(ans_abc, main = "No priors ABC")
#' plot(ans_nr_dbeta, main = "L-BFGS-B w/ Prior for Psi ~ beta(2,10)")
#' plot(ans_abc_dbeta, main = "ABC w/ Prior for Psi ~ beta(2,10)")
#' 
#' par(oldpar)
#' 
#' 
#' # Using the MCMC ------------------------------------------------------------
#' 
#' \dontrun{
#' 
#' set.seed(1233)
#' # Simulating a tree
#' tree <- sim_tree(200)$edges
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
#' ans_mcmc <- phylo_mcmc(
#'   rep(.1, 5), dat,
#'   control = list(nbatch = 2e5, burnin=1000, thin=200, scale=2e-2)
#' )
#' }
#' 
#' @name mle
NULL

#' @rdname mle
#' @export
phylo_mle <- function(
  dat,
  method        = "L-BFGS-B",
  priors        = NULL, 
  control       = list(),
  params        = rep(.05, 5),
  lower         = 0,
  upper         = 1,
  fix.params    = c(psi0 = FALSE, psi1 = FALSE, mu0 = FALSE, mu1 = FALSE, Pi = FALSE)
) {
  
  
  # Checking params
  if (length(params) != 5)
    stop("-params- must be of length 5.")
  
  if (!is.numeric(params))
    stop("-params- must be a numeric vector")
  
  # In case of fixing parameters
  par0 <- params
  
  # Both available for ABC
  control$fnscale  <- -1
  control$parscale <- rep(1, 5)
  
  # Creating the objective function
  fun <- if (length(priors)) {
    function(params, dat) {
      
      # Checking whether params are fixed or not
      params <- ifelse(fix.params, par0, params)
      
      psi <- params[1:2] 
      mu  <- params[3:4] 
      Pi  <- params[5]
      
      ll <- LogLike(
        annotations = dat$annotations, 
        offspring   = dat$offspring,
        noffspring  = dat$noffspring, 
        psi         = psi, 
        mu          = mu, 
        Pi          = Pi, 
        verb_ans    = FALSE, 
        check_dims  = FALSE
        )$ll + sum(log(priors(params)))
      
      # Checking if we got a finite result
      if (is.infinite(ll)) return(.Machine$double.xmax*sign(ll)*1e-10)
      ll
    }
  } else {
    function(params, dat) {
      # Checking whether params are fixed or not
      params <- ifelse(fix.params, par0, params)
      
      psi <- params[1:2] 
      mu  <- params[3:4] 
      Pi  <- params[5]
      
      ll <- LogLike(
        annotations = dat$annotations, 
        offspring   = dat$offspring,
        noffspring  = dat$noffspring, 
        psi         = psi, 
        mu          = mu, 
        Pi          = Pi, 
        verb_ans    = FALSE, 
        check_dims  = FALSE
      )$ll

      # Checking if we got a finite result
      if (is.infinite(ll)) return(.Machine$double.xmax*sign(ll)*1e-10)
      ll
    }
  }
  
  # Optimizing
  if (method == "ABC") {
    # Checking ABC args
    if (!length(control$lb))       control$lb       <- lower
    if (!length(control$ub))       control$ub       <- upper
    if (!length(control$maxCycle)) control$maxCycle <- 500L
    if (!length(control$criter))   control$criter   <- 50L

    ans <-
      do.call(ABCoptim::abc_cpp, c(list(par = params, fun, dat=dat), control))
    ans$convergence <- NA
    ans$message     <- NA
  } else {
    
    # # Defining gradient
    # grd <- function(params, dat) {
    #   numDeriv::grad(fun, params, dat=dat, method.args = list(v = 6))
    # }

    # Try to solve it
    ans <- do.call(
      stats::optim, 
      c(
        list(par = params, fn = fun, method=method, upper = upper, lower = lower,
             hessian=FALSE, dat = dat, control=control)
        )
      )
    
    ans     <- list(
      hist  = rbind(params, ans$par),
      par   = ans$par,
      value = ans$value,
      convergence = ans$convergence,
      message = ans$message,
      counts = ans$counts["function"]
      )
    
    
  }
  
  # Computing the hessian (information matrix)
  hessian <- stats::optimHess(ans$par, fun, dat = dat, control = control)
  
  # Working on answer
  env <- new.env()
  environment(fun)    <- env
  environment(dat)    <- env
  environment(par0)   <- env
  environment(fix.params) <- env
  if (length(priors)) 
    environment(priors) <- env
  
  # Naming the parameter estimates
  names(ans$par) <- c("psi0", "psi1", "mu0", "mu1", "Pi")
  
  # Hessian for observed information matrix
  dimnames(hessian) <- list(names(ans$par), names(ans$par))
  
  # Returning
  structure(
    list(
      par        = ans$par,
      hist       = ans$hist,
      ll         = ans$value,
      counts     = ans$counts,
      convergence = ans$convergence,
      message    = ans$message,
      fun        = fun,
      priors     = priors,
      dat        = dat,
      par0       = par0,
      fix.params = fix.params,
      method     = method,
      varcovar   = solve(-hessian, tol = 1e-100)
    ),
    class = "phylo_mle"
  )
}

#' @export
#' @rdname mle
print.phylo_mle <- function(x, ...) {
  # Function to print a bar with variable width
  catbar <- function() paste0(rep("-",options()$width), collapse="")
  
  sderrors   <- sqrt(diag(x$varcovar))
  props      <- with(x$dat, table(annotations[noffspring == 0,,drop=FALSE]))
  propspcent <- prop.table(props)*100
  
  with(x, {
    cat(
      sep = "\n",
      "ESTIMATION OF ANNOTATED PHYLOGENETIC TREE",
      sprintf(
        "ll: %9.4f,\nMethod used: %s (%i iterations)", ll, method, x$counts),
      if (method %in% c("mcmc", "ABC")) 
        NULL
      else
        sprintf("convergence: %i (see ?optim)", convergence)
      ,
      sprintf("Leafs\n # of Functions %i", ncol(dat$annotations)),
      paste0(sprintf(" # of %s: %5i (%2.0f%%)", names(props), props, propspcent), collapse="\n"),
            "\n         Estimate  Std. Error",
      sprintf(" psi[0]    %6.4f      %6.4f", par["psi0"], sderrors["psi0"]),
      sprintf(" psi[1]    %6.4f      %6.4f", par["psi1"], sderrors["psi1"]),
      sprintf(" mu[0]     %6.4f      %6.4f", par["mu0"], sderrors["mu0"]),
      sprintf(" mu[1]     %6.4f      %6.4f", par["mu1"], sderrors["mu1"]),
      sprintf(" Pi        %6.4f      %6.4f", par["Pi"], sderrors["Pi"])
      )
    
  })
  
  invisible(x)
}

#' @export
#' @rdname mle
coef.phylo_mle <- function(object, ...) {
  object$par
}

#' @export
#' @rdname mle
vcov.phylo_mle <- function(object, ...) {
  object$varcovar
}

#' @param x An object of class \code{phylo_mle}.
#' @param y Ignored.
#' @param main Passed to plot.
#' @param xlab Passed to plot.
#' @param ylab Passed to plot.
#' @param type Passed to plot.
#' @param ... Further arguments passed to the method.
#' @param addlegend Logical scalar. When \code{TRUE} adds extra info to the
#' plot as a legend including ll at the optimum and parameter values.
#' @rdname mle
#' @export
plot.phylo_mle <- function(
  x,
  y = NULL,
  main = "",
  xlab = "Step",
  ylab = "Log-Likelihood",
  type = ifelse(x$method %in% c("mcmc", "ABC"), "l", "p"),
  addlegend = TRUE,
  ...
  ) {
  
  # Filtering complete cases
  # x$hist <- x$hist[apply(x$hist, 1, function(x) !any(is.na(x))),]
  if (x$method == "mcmc")
    x$hist <- do.call(rbind, x$hist)
  
  plot(
    y = apply(as.matrix(x$hist), 1, x$fun, dat=x$dat),
    x = 1:nrow(x$hist),
    type = type,
    main = main,
    xlab = xlab,
    ylab = ylab,
    ...
    )

  if (addlegend) {
    
    numbers <- sprintf("%.5f", with(x, c(ll, par)))
    numbers <- c(
      bquote(L(theta~"|"~X) == .(numbers[1])),
      bquote(psi[0] == .(numbers[2])),
      bquote(psi[1] == .(numbers[3])),
      bquote(mu[0] == .(numbers[4])),
      bquote(mu[1] == .(numbers[5])),
      bquote(pi == .(numbers[6]))
    )
    
    legend("bottomright",legend = sapply(numbers, as.expression),bty = "n")
  }
  
  
}

#' @rdname mle
#' @return In the case of \code{phylo_mcmc}, \code{hist} is an object of class
#' \code{\link[coda:mcmc.list]{mcmc.list}}.
#' @export
phylo_mcmc <- function(
  params,
  dat,
  priors        = NULL,
  fix.params    = c(psi0 = FALSE, psi1 = FALSE, mu0 = FALSE, mu1 = FALSE, Pi = FALSE),
  control       = list()
) {
  
  # Checking control
  if (!length(control$nbatch)) control$nbatch <- 2e3
  if (!length(control$scale))  control$scale  <- .01
  if (!length(control$ub))     control$ub     <- rep(1, 5)
  if (!length(control$lb))     control$lb     <- rep(0, 5)
  
  # Checking params
  if (length(params) != 5)
    stop("-params- must be of length 5.")
  
  if (!is.numeric(params))
    stop("-params- must be a numeric vector")
  
  # In case of fixing parameters
  par0 <- params
  
  # Creating the objective function
  fun <- if (length(priors)) {
    function(params, dat) {

      # Checking whether params are fixed or not
      params <- ifelse(fix.params, par0, params)
      
      LogLike(
        annotations = dat$annotations,
        offspring   = dat$offspring,
        noffspring  = dat$noffspring,
        psi         = params[1:2] ,
        mu          = params[3:4] ,
        Pi          = params[5],
        verb_ans    = FALSE, 
        check_dims  = FALSE
      )$ll + sum(log(priors(params)))
      
    }
  } else {
    function(params, dat) {

      # Checking whether params are fixed or not
      params <- ifelse(fix.params, par0, params)
      
      LogLike(
        annotations = dat$annotations,
        offspring   = dat$offspring,
        noffspring  = dat$noffspring,
        psi         = params[1:2] ,
        mu          = params[3:4] ,
        Pi          = params[5],
        verb_ans    = FALSE, 
        check_dims  = FALSE
      )$ll
      
    }
  }
  
  # Naming the parameter estimates
  names(params) <- c("psi0", "psi1", "mu0", "mu1", "Pi")
  
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
  environment(fix.params) <- env
  if (length(priors)) 
    environment(priors) <- env
  
  # Returning
  structure(
    list(
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
      fix.params = fix.params,
      method     = "mcmc",
      varcovar   = var(do.call(rbind, ans))
    ),
    class = "phylo_mle"
  )
}

