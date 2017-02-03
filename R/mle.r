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
#' \code{maxiter} \tab item Integer scalar. Maximum number of steps in the Newton-Raphson
#' algorithm. Default \code{20L}.\cr
#' \code{criter} \tab Numeric scalar. Stoping criteria for the Newton-Raphson
#' algorithm. Default \code{1e-15}.\cr
#' \code{method.args} \tab A list of arguments passed to
#' \code{\link[numDeriv:jacobian]{hessian,jacobian}} from the \CRANpkg{numDeriv}
#' package. Default \code{list(d = .0001)}\cr
#' \code{solve.tol} \tab Numeric scalar passed to \code{\link{solve}}. Default
#' \code{1e-40}.
#' }
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
#' # Using the package data ----------------------------------------------------
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
#' # Computing Estimating the parameters 
#' ans_nr  <- phylo_mle(rep(.5,5), O)
#' ans_abc <- phylo_mle(rep(.5,5), O, useABC = TRUE)
#' 
#' # Plotting the path
#' with(ans_nr, plot(
#'   - apply(hist, 1, fun),
#'   type = "l",
#'   xlab = "Step",
#'   ylab = "Log-Likelihood"
#' ))
#' 
#' # Computing Estimating the parameters Using Priors for PSI 
#' mypriors <- function(params) {
#'     dbeta(params[1:2], 2, 10)
#' }
#' ans_nr_dbeta <- phylo_mle(rep(.5,5), O, priors = mypriors)
#' ans_abc_dbeta <- phylo_mle(rep(.5,5), O, priors = mypriors, useABC = TRUE)
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
#' }
#' # Adding some zeros to the data ---------------------------------------------
#' set.seed(1231)
#' mysample <- which(rowSums(experiment[,1:3]) == 27)
#' mysample <- sample(mysample, 20)
#' experiment[mysample,1:3] <- 0
#' 
#' O <- get_offspring(
#'   experiment, "LeafId", 
#'   tree, "NodeId", "ParentId"
#' )
#' 
#' ans_nr  <- phylo_mle(rep(.5,5), O)
#' ans_abc <- phylo_mle(rep(.5,5), O, useABC = TRUE)
#' 
#' 
#' # Computing Estimating the parameters Using Priors for PSI ------------------
#' mypriors <- function(params) {
#'   dbeta(params[1:2], 2, 10)
#' }
#' ans_nr_dbeta <- phylo_mle(rep(.5,5), O, priors = mypriors)
#' ans_abc_dbeta <- phylo_mle(rep(.5,5), O, priors = mypriors, useABC = TRUE)
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
  
  # Auxiliary functions for 
  expit <- function(x) exp(x)/(1 + exp(x))
  
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
    
    if (!length(control$maxiter))     control$maxiter    <- 20L
    if (!length(control$criter))      control$criter     <- 1e-15
    if (!length(control$solve.tol))   control$solve.tol  <- 1e-40
    if (!length(control$method.args)) control$method.args  <- list(d = .0001)
    
    # Creating space
    PARAMS <- matrix(ncol = 5, nrow = control$maxiter)
    
    # Newton-Raphson. Observe that in each evaluation of -fun- we apply the
    # -expit- function to the parameters so that those are transformed to [0,1]
    for (i in 1:control$maxiter) {
      
      PARAMS[i,] <- params
      
      # Computing jacobian and hessian
      fun_jacb   <- numDeriv::jacobian(fun, expit(params), method.args = control$method.args)
      fun_hess   <- numDeriv::hessian(fun, expit(params), method.args = control$method.args)
      
      # Updating step
      params <- params - fun_jacb %*% solve(fun_hess, tol = control$solve.tol) 
      
      # Error
      if (is.na(fun(expit(params)))) {
        message("Undefined value of fun(params).")
        break
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
  
  # Returning
  structure(
    list(
      par = ans$par,
      hist = ans$hist,
      value = ans$value,
      ll    = -ans$value,
      fun = fun,
      priors = priors,
      dat = dat,
      par0 = par0,
      fix.params = fix.params,
      method = ifelse(useABC, "ABC", "NR")
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
  
  # Filtering complete cases
  x$hist <- x$hist[apply(x$hist, 1, function(x) !any(is.na(x))),]
  
  with(x,
       plot(
         ifelse(method == "mcmc", 1, -1)*apply(hist, 1, fun),
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
