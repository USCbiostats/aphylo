#' Maximum Likelihood Estimation of the Phylogenetic Tree
#'
#' The optimization is done either via \code{abc_cpp} or Newton-Raphson.
#'
#' @param params A vector of length 5 with initial parameters. In particular
#' \code{psi[1]}, \code{psi[2]}, \code{mu[1]}, \code{mu[2]}, and \code{Pi}.
#' @param dat An object of class \code{phylo_offspring} as returned by
#' \code{\link{get_offspring}}.
#' @param useABC Logical scalar. When \code{TRUE}, uses Artificial Bee Colony
#' optimization algorithm instead. 
#' @param priors A list of length 3 with functions named \code{psi}, \code{mu},
#' \code{Pi}
#' @param fix.params A Logical vector of length 5. Whether or not to fix
#' a particular parameter to that of what was specified in \code{params}.
#' @param control A list with parameters for the optimization method (see
#' details).
#' @param object An object of class \code{phylo_mle}.
#' 
#' @details When \code{useABC = TRUE}, the optimization is done via 
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
#' The default \code{control} parameters for \bold{Newton-Raphson method} are
#' the following:
#' 
#' \tabular{ll}{
#' \code{maxiter} \tab Integer scalar. Maximum number of steps in the Newton-Raphson
#' algorithm. Default \code{100L}.\cr
#' \code{criter} \tab Numeric scalar. Stoping criteria for the Newton-Raphson
#' algorithm. Default \code{1e-15}.\cr
#' \code{method.args} \tab A list of arguments passed to
#' \code{\link[numDeriv:jacobian]{hessian,jacobian}} from the \CRANpkg{numDeriv}
#' package. Default \code{list(d = .0001)}\cr
#' \code{solve.tol} \tab Numeric scalar passed to \code{\link{solve}}. Default
#' \code{1e-40}.
#' }
#' 
#' Notice that the algorithm is somewhat numerically unstable when using starting points
#' with values higher than \code{0.4}. It is recomended to use values closer to
#' \code{0.1} instead.
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
#' \item{hist}{A numeric matrix of size \code{counts*5} with the solution path.}
#' \item{value}{A numeric scalar with the value of \code{fun(par)}}
#' \item{ll}{A numeric scalar with the value of \code{-fun(par)}}
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
#' dat <- sim_annotated_tree(250, P=2)
#' 
#' 
#' # Computing Estimating the parameters 
#' # for some reason, starting with parameters equal to .5 breaks NR.
#' ans_nr  <- phylo_mle(rep(.1,5), dat)
#' ans_abc <- phylo_mle(rep(.1,5), dat, useABC = TRUE)
#' 
#' # Plotting the path
#' with(ans_nr, plot(
#'  - apply(hist, 1, fun),
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
#' ans_nr_dbeta <- phylo_mle(rep(.1,5), dat, priors = mypriors)
#' ans_abc_dbeta <- phylo_mle(rep(.1,5), dat, priors = mypriors, useABC = TRUE)
#' 
#' # Plotting the path
#' oldpar <- par(no.readonly = TRUE)
#' par(mfrow = c(2, 2))
#' 
#' plot(ans_nr, main = "No priors NR")
#' plot(ans_abc, main = "No priors ABC")
#' plot(ans_nr_dbeta, main = "NR w/ Prior for Psi ~ beta(2,10)")
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
#' ans_mcmc <- phylo_mcmc(
#'   rep(.5, 5), dat,
#'   control = list(nbatch = 2e5, burnin=1000, thin=200, scale=2e-2)
#' )
#' }
#' 
#' @name mle
NULL

#' @rdname mle
#' @export
phylo_mle <- function(
  params,
  dat,
  useABC        = FALSE,
  priors        = NULL, 
  control       = list(),
  fix.params    = c(psi0 = FALSE, psi1 = FALSE, mu0 = FALSE, mu1 = FALSE, Pi = FALSE)
) {
  
  
  # Checking params
  if (length(params) != 5)
    stop("-params- must be of length 5.")
  
  if (!is.numeric(params))
    stop("-params- must be a numeric vector")
  
  
  
  # In case of fixing parameters
  par0 <- params
  
  # Creating the objective function
  fun <- if (length(priors)) {
    function(params) {
      
      # Checking whether params are fixed or not
      params <- ifelse(fix.params, par0, params)
      
      psi <- params[1:2] 
      mu  <- params[3:4] 
      Pi  <- params[5]
      Pi  <- c(1 - Pi, Pi)
      
      - LogLike(dat$experiment, dat$offspring, dat$noffspring, psi, mu, Pi, FALSE)$ll -
         sum(log(priors(params)))

    }
  } else {
    function(params) {
      # Checking whether params are fixed or not
      params <- ifelse(fix.params, par0, params)
      
      psi <- params[1:2] 
      mu  <- params[3:4] 
      Pi  <- params[5]
      Pi  <- c(1 - Pi, Pi)
      
      - LogLike(dat$experiment, dat$offspring, dat$noffspring, psi, mu, Pi, FALSE)$ll
    }
  }
  
  # Optimizing
  if (useABC) {
    # Checking ABC args
    if (!length(control$lb))       control$lb       <- 1e-20
    if (!length(control$ub))       control$ub       <- 1 - 1e-20
    if (!length(control$maxCycle)) control$maxCycle <- 500L
    if (!length(control$criter))   control$criter   <- 50L

    ans <-
      do.call(ABCoptim::abc_cpp, c(list(par = params, fn = fun), control))
  } else {
    
    # Checking NR args
    if (!length(control$maxiter))     control$maxiter    <- 100L
    if (!length(control$criter))      control$criter     <- 1e-15
    if (!length(control$solve.tol))   control$solve.tol  <- 1e-40
    if (!length(control$method.args)) control$method.args  <- list(d = .0001)
    
    # Auxiliary functions for 
    expit <- function(x) exp(x)/(1 + exp(x))
    logit <- function(x) - log(1/x - 1)
    
    # Creating space
    PARAMS     <- matrix(ncol = 5, nrow = control$maxiter)
    params     <- logit(params)
    
    # Newton-Raphson. Observe that in each evaluation of -fun- we apply the
    # -expit- function to the parameters so that those are transformed to [0,1]
    for (i in 1:control$maxiter) {
      
      # Current value
      PARAMS[i,] <- params
      
      # Computing jacobian and hessian
      fun_jacb <- numDeriv::jacobian(fun, expit(params), method.args = control$method.args)
      fun_hess <- numDeriv::hessian(fun, expit(params), method.args = control$method.args)
      
      # Updating step
      params <- params - fun_jacb %*% solve(fun_hess, tol = control$solve.tol) 

      # Error
      if (is.na(fun(expit(params)))) {
        stop("Undefined value of fun(params). Try using a different starting point for -params-.")
      }
        
      # Stopping criteria
      val <- abs(fun(expit(params)) - fun(expit(PARAMS[i, ])))
      if (val < control$criter) {
        break
      }
    }
    
    # Wrapping value
    ans <- list(
      hist  = expit(PARAMS[1:i,,drop = FALSE]),
      par   = expit(as.vector(params)),
      value = fun(expit(params))
    )
  }
  
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
  hessian <- numDeriv::hessian(function(p) -fun(p), ans$par, method.args = control$method.args)
  dimnames(hessian) <- list(
    names(ans$par),
    names(ans$par)
  )
  
  # Returning
  structure(
    list(
      par        = ans$par,
      hist       = ans$hist,
      value      = ans$value,
      ll         = -ans$value,
      fun        = fun,
      priors     = priors,
      dat        = dat,
      par0       = par0,
      fix.params = fix.params,
      method     = ifelse(useABC, "ABC", "NR"),
      varcovar   = -solve(hessian)
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
  props      <- with(x$dat, table(experiment[noffspring == 0,,drop=FALSE]))
  propspcent <- prop.table(props)*100
  
  with(x, {
    cat(
      sep = "\n",
      "ESTIMATION OF ANNOTATED PHYLOGENETIC TREE",
      sprintf(
        "ll: %9.4f,\nMethod used: %s (%i steps)\nLeafs\n # of Functions %i",
        ll, method, nrow(x$hist), ncol(x$dat$experiment)
        ),
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
  type = "l",
  addlegend = TRUE,
  ...
  ) {
  
  # Filtering complete cases
  x$hist <- x$hist[apply(x$hist, 1, function(x) !any(is.na(x))),]
  
  with(x,
       plot(
         ifelse(method == "mcmc", 1, -1)*apply(as.matrix(hist), 1, fun),
         type = type,
         main = main,
         xlab = xlab,
         ylab = ylab,
         ...
       ))

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
