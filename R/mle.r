#' Maximum Likelihood Estimation of the Phylogenetic Tree
#'
#' The optimization is done either via \code{abc_cpp} or Newton-Raphson.
#'
#' @param params A vector of length 5 with initial parameters. In particular
#' \code{psi[1]}, \code{psi[2]}, \code{mu[1]}, \code{mu[2]}, and \code{Pi}.
#' @param dat An object of class \code{phylo_offspring} as returned by
#' \code{\link{get_offspring}}.
#' @param maxiter Integer scalar. Maximum number of steps in the Newton-Raphson
#' algorithm.
#' @param criter Numeric scalar. Stoping criteria for the Newton-Raphson algorithm.
#' @param useABC Logical scalar. When \code{TRUE}, uses Artificial Bee Colony
#' optimization algorithm instead. 
#' @param priors A list of length 3 with functions named \code{psi}, \code{mu},
#' \code{Pi}
#' @param abcoptim.args A list of arguments to be passed to
#' \code{\link[ABCoptim:abc_cpp]{abc_cpp}} in the \pkg{ABCoptim} package.
#' 
#' @return 
#' A list of class \code{phylo_mle} with the following elements:
#' \item{par}{A numeric vector of length 5 with the solution.}
#' \item{hist}{A numeric matrix of size \code{counts*5} with the solution path.}
#' \item{value}{A numeric scalar with the value of \code{fun(par)}}
#' \item{fun}{A function (the objective function).}
#' \item{dat}{The data \code{dat} provided to the function.}
#' 
#' @examples 
#' 
#' \dontrun{
#' # Loading data
#' data(experiment)
#' data(tree)
#' 
#' # Preprocessing the data
#' O <- get_offspring(
#'   experiment, "LeafId", 
#'   tree, "NodeId", "ParentId"
#' )
#' 
#' # Computing Estimating the parameters ---------------------------------------
#' ans_nr  <- mle(rep(.5,5), O)
#' ans_abc <- mle(rep(.5,5), O, useABC = TRUE)
#' 
#' # Plotting the path
#' with(ans_nr, plot(
#'   - apply(hist, 1, fun),
#'   type = "l",
#'   xlab = "Step",
#'   ylab = "Log-Likelihood"
#' ))
#' 
#' # Computing Estimating the parameters Using Priors for PSI ------------------
#' mypriors <- function(params) {
#'     dbeta(params[1:2], 1, 9)
#' }
#' ans_nr_dbeta <- mle(rep(.5,5), O, priors = mypriors)
#' 
#' # Plotting the path
#' oldpar <- par(no.readonly = TRUE)
#' par(mfrow = c(2, 2))
#' 
#' plot(ans_nr, main = "No priors NR")
#' plot(ans_abc, main = "No priors ABC")
#' plot(ans_nr_dbeta, main = "Prior for Psi ~ beta(1,9)")
#' 
#' par(oldpar)
#' }
#' @name mle
NULL

#' @rdname mle
#' @export
mle <- function(
  params,
  dat,
  maxiter       = 20L,
  criter        = 1e-5,
  useABC        = FALSE,
  priors        = NULL, 
  abcoptim.args = list(ub = 1, lb = 0, maxCycle = 500L, criter = 50L)
) {
  
  # Checking params
  if (length(params) != 5)
    stop("-params- must be of length 5.")
  
  if (!is.numeric(params))
    stop("-params- must be a numeric vector")
  
  # Auxiliary functions for 
  expit <- function(x) exp(x)/(1 + exp(x))
  logit <- function(x) log(x/(1 - x))
  
  # Creating the objective function
  fun <- if (length(priors)) {
    function(params) {
      psi <- params[1:2] 
      mu  <- params[3:4] 
      Pi  <- params[5]
      Pi  <- c(1 - Pi, Pi)
      
      - LogLike(dat$experiment, dat$offspring, dat$noffspring, psi, mu, Pi, FALSE)$ll - 
        sum(log(priors(params)))
    }
  } else {
    function(params) {
      psi <- params[1:2] 
      mu  <- params[3:4] 
      Pi  <- params[5]
      Pi  <- c(1 - Pi, Pi)
      
      - LogLike(dat$experiment, dat$offspring, dat$noffspring, psi, mu, Pi, FALSE)$ll
    }
  }
  
  # Optimizing
  if (useABC) {
    ans <-
      do.call(ABCoptim::abc_cpp, c(list(par = params, fn = fun), abcoptim.args))
  } else {
    
    # Creating space
    PARAMS <- matrix(ncol = 5, nrow = maxiter)
    
    # Newton-Raphson. Observe that in each evaluation of -fun- we apply the
    # -expit- function to the parameters so that those are transformed to [0,1]
    for (i in 1:maxiter) {
      
      PARAMS[i,] <- params
      
      # Computing jacobian and hessian
      fun_jacb   <- numDeriv::jacobian(fun, expit(params), method.args = list(d = .025))
      fun_hess   <- numDeriv::hessian(fun, expit(params), method.args = list(d = .025))
      
      # Updating step
      params <- params - fun_jacb %*% solve(fun_hess, tol = 1e-40)
      
      # Error
      if (is.na(fun(expit(params)))) {
        message("Undefined value of fun(params).")
        break
      }
        
      # Stopping criteria
      val <- abs(fun(expit(params)) - fun(expit(PARAMS[i, ])))
      if (val < criter) {
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
  if (length(priors)) 
    environment(priors) <- env
  
  # Naming the parameter estimates
  names(ans$par) <- c("psi0", "psi1", "mu0", "mu1", "Pi")
  
  # Returning
  structure(
    list(
      par = ans$par,
      hist = ans$hist,
      value = ans$value,
      ll    = -ans$value,
      fun = fun,
      priors = priors,
      dat = dat
    ),
    class = "phylo_mle"
  )
}

#' @param x An object of class \code{phylo_mle}.
#' @param y Ignored.
#' @param main Passed to plot.
#' @param xlab Passed to plot.
#' @param ylab Passed to plot.
#' @param type Passed to plot.
#' @param ... Further arguments passed to plot
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
  
    with(x,
         plot(
           -apply(hist, 1, fun),
           type = type,
           main = main,
           xlab = xlab,
           ylab = ylab,
           ...
         ))
  
  if (addlegend) {
    
    numbers <- sprintf("%.4f", with(x, c(ll, par)))
    numbers <- c(
      bquote(L(theta || X) == .(numbers[1])),
      bquote(psi[0] == .(numbers[2])),
      bquote(psi[1] == .(numbers[3])),
      bquote(mu[0] == .(numbers[4])),
      bquote(mu[1] == .(numbers[5])),
      bquote(pi == .(numbers[6]))
    )
    
    legend("bottomright",legend = sapply(numbers, as.expression),bty = "n")
  }
  
  
  }
