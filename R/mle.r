#' Maximum Likelihood Estimation of the Phylogenetic Tree.
#' @param params A vector of length 5 with initial parameters. In particular
#' \code{psi[1]}, \code{psi[2]}, \code{mu[1]}, \code{mu[2]}, and \code{Pi}.
#' @param dat An object of class \code{phylo_offspring}.
#' @param prior A list of length 3 with functions named \code{psi}, \code{mu},
#' \code{Pi}
#' @param abcoptim.args A list of arguments to be passed to
#' \code{\link[ABCoptim:abc_cpp]{abc_cpp}} in the \pkg{ABCoptim} package.
#' @examples 
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
#' ans0 <- mle(rep(.5,5), O)
#' 
#' # Plotting the path
#' with(ans0, plot(
#'   - apply(abc$hist, 1, fun),
#'   type = "l",
#'   xlab = "Step",
#'   ylab = "Log-Likelihood"
#' ))
#' 
#' # Computing Estimating the parameters Using Priors for PSI ------------------
#' ans1 <- mle(rep(.5,5), O,
#'     priors        = list(psi = function(x) dbeta(x, 1, 9)))
#' 
#' # Plotting the path
#' oldpar <- par(no.readonly = TRUE)
#' par(mfrow = c(1, 2))
#' 
#' plot(ans0, main = "No priors")
#' plot(ans0, main = "Prior for Psi ~ beta(1,9)")
#' 
#' par(oldpar)
#' @export
mle <- function(
  params,
  dat,
  maxiter       = 20L,
  criter        = 1e-5,
  useABC        = FALSE,
  priors        = list(psi = NULL, mu = NULL, Pi = NULL), 
  abcoptim.args = list(ub = .9999, lb = .0001, maxCycle = 500L, criter = 50L)
) {
  
  # Checking params
  if (length(params) != 5)
    stop("-params- must be of length 5.")
  
  if (!is.numeric(params))
    stop("-params- must be a numeric vector")
  
  # Checking priors
  priormult <- "function(psi, mu, Pi) { 0"
  for (p in c("psi", "mu", "Pi")) {
    # If there is no prior
    if (length(priors[[p]]) == 0) {
      priors[[p]] <- function(x) x
      
    } else {
      # If there's any prior
      priormult <- paste(
        priormult,
        " + log(priors[['",p,"']](", p, "[1]))",
        " + log(priors[['",p,"']](", p, "[2]))",
        sep = "")
    }
  }
  
  # Coercing into a function
  priormult <- paste(priormult, "}")
  priormult <- eval(parse(text = priormult))
  
  # Auxiliary functions for 
  expit <- function(x) exp(x)/(1 + exp(x))
  logit <- function(x) log(x/(1 - x))
  
  # Optimizing
  if (useABC) {
    fun <- function(params) {
      psi <- params[1:2] # c(0.020,0.010)
      mu  <- params[3:4] # c(0.004,.001)
      Pi  <- params[5]
      Pi  <- c(1 - Pi, Pi) # c(1-0.1,.1) 
      
      - LogLike(dat$experiment, dat$offspring, dat$noffspring, psi, mu, Pi)$ll - 
        priormult(psi, mu, Pi)
    }
    
    ans <- do.call(ABCoptim::abc_cpp,
                   c(list(par = params, fn = fun), abcoptim.args))
  } else {
    
    # Objective function
    fun <- function(params) {
      params <- expit(params)
      
      psi <- params[1:2] # c(0.020,0.010)
      mu  <- params[3:4] # c(0.004,.001)
      Pi  <- params[5]
      Pi  <- c(1 - Pi, Pi) # c(1-0.1,.1) 
      
      - LogLike(dat$experiment, dat$offspring, dat$noffspring, psi, mu, Pi)$ll - 
        priormult(psi, mu, Pi)
    }
    
    # Creating space
    PARAMS <- matrix(ncol = 5, nrow = maxiter)
    params <- logit(params)
    
    # Looping
    for (i in 1:maxiter) {
      
      PARAMS[i,] <- params
      params0    <- params
      
      # Computing jacobian and hessian
      fun_jacb   <- numDeriv::jacobian(fun, params, method.args = list(d = .025))
      fun_hess   <- numDeriv::hessian(fun, params, method.args = list(d = .025))
      
      # Updating step
      params <- params - fun_jacb %*% solve(fun_hess, tol = 1e-40)
      
      # Error
      if (is.na(fun(params))) {
        message("Undefined value of fun(params).")
        break
      }
        
      
      # Stopping criteria
      val <- abs(fun(params) - fun(PARAMS[i, ]))
      if (val < criter) {
        break
      }
    }
    
    # Wrapping value
    ans <- list(
      hist  = expit(PARAMS[1:i,,drop = FALSE]),
      par   = expit(params),
      value = -fun(expit(params))
    )
  }
  
  # Working on answer
  env <- new.env()
  environment(fun)       <- env
  environment(priormult) <- env
  environment(dat)       <- env
  
  # Returning
  structure(list(
    abc = ans,
    priormult = priormult,
    fun = fun,
    dat = dat
  ), class = "phylo_mle")
}

#' @param x An object of class \code{phylo_mle}.
#' @param y Ignored.
#' @param main Passed to plot.
#' @param xlab Passed to plot.
#' @param ylab Passed to plot.
#' @param type Passed to plot.
#' @param ... Further arguments passed to plot
#' @rdname mle
#' @export
plot.phylo_mle <- function(
  x,
  y = NULL,
  main = "Prior for Psi ~ beta(1,9)",
  xlab = "Step",
  ylab = "Log-Likelihood",
  type = "l",
  ...
  ) {
  
    with(x,
         plot(
           -apply(abc$hist, 1, fun),
           type = type,
           main = main,
           xlab = xlab,
           ylab = ylab,
           ...
         ))
  }
