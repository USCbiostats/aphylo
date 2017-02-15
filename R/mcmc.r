#' Markov Chain Monte Carlo
#' @param fun A function. Returns the log-likelihood
#' @param initial A numeric vector. Initial values of the parameters.
#' @param nbatch Integer scalar. Number of MCMC runs.
#' @param thin Integer scalar. Passed to \code{\link[coda:mcmc]{coda::mcmc}}.
#' @param scale Numeric scalar. Step size, \eqn{z\times scale}{z*scale}, where
#' \eqn{z\sim N(0,1)}{z~N(0,1)}.
#' @param burnin Integer scalar. Number of burn-in samples. Passed to 
#' \code{\link[coda:mcmc]{coda::mcmc}} as \code{init}.
#' @param lb Numeric vector of length \code{length(initial)}. Lower bounds
#' @param ub Numeric vector of length \code{length(initial)}. Upper bounds
#' @param inf.eps Numeric scalar. See details.
#' 
#' @details This function implements MCMC using the Metropolis-Hastings ratio with
#' scaled standard normal propositions for each parameter. For each parameter
#' the transition function is 
#' 
#' \deqn{
#' \theta' = \theta + scale*z
#' }
#' 
#' Where \eqn{z} has standard normal distribution. The MCMC follows a block
#' sampling scheme, i.e. proposed states are either accepted or rejected
#' altogether.
#' 
#' Lower and upper bounds are treated using reflecting boundaries, this is, 
#' if the proposed \eqn{\theta'} is greater than the \code{ub}, then \eqn{\theta' - ub}
#' is substracted from \eqn{ub}. At the same time, if it is less than \code{lb}, then
#' \eqn{lb - \theta'} is added to \code{lb} iterating until \eqn{\theta} is within
#' \code{[lb, ub]}.
#' 
#' If \code{name(initial) == NULL}, then a names in the form of \code{par1, par2, ...}
#' will be assigned to the variables.
#' 
#' If \code{fun} returns \code{-Inf} or \code{Inf}, instead of discarding such parameters
#' or accept/reject with certainty, the function replaces the values to
#' \code{-.Machine$double.xmax*inf.eps} and \code{.Machine$double.xmax*inf.eps}
#' respectively so that the algorithm is aperiodic.
#' 
#' @return An object of class \code{\link[coda:mcmc]{mcmc}} from the \CRANpkg{coda}
#' package. The \code{mcmc} object is a matrix with one column per parameter,
#' and \code{nbatch} rows.
#' 
#' 
#' @export
#' @examples 
#' # Parameters
#' set.seed(1231)
#' n <- 1e3
#' pars <- c(mean = 2.6, sd = 3)
#' 
#' # Generating data and writing the log likelihood function
#' D <- rnorm(n, pars[1], pars[2])
#' fun <- function(x) {
#'   x <- log(dnorm(D, x[1], x[2]))
#'   if (any(is.infinite(x)))
#'     return(-Inf)
#'   sum(x)
#' }
#' ans <- MCMC(fun, c(mu=1, sigma=1), nbatch = 2e3, scale = .1, ub = 10, lb = 0)
#' oldpar <- par(no.readonly = TRUE)
#' par(mfrow = c(1,2))
#' boxplot(as.matrix(ans), 
#'         main = expression("Posterior distribution of"~mu~and~sigma),
#'         names =  expression(mu, sigma), horizontal = TRUE,
#'         col  = blues9[c(4,9)],
#'         sub = bquote(mu == .(pars[1])~", and"~sigma == .(pars[2]))
#' )
#' abline(v = pars, col  = blues9[c(4,9)], lwd = 2, lty = 2)
#' 
#' plot(apply(as.matrix(ans), 1, fun), type = "l",
#'      main = "LogLikelihood",
#'      ylab = expression(L("{"~mu,sigma~"}"~"|"~D)) 
#' )
#' par(oldpar)
#' @aliases Metropolis-Hastings
MCMC <- function(
  fun,
  initial, 
  nbatch = 2e4L,
  thin   = 1L,
  scale  = 1,
  burnin = 1e3L,
  ub = rep(.Machine$double.xmax, length(initial)),
  lb = rep(-.Machine$double.xmax, length(initial)),
  inf.eps = 1e-100
  ) {
  
  # Adding names
  cnames <- names(initial)
  if (!length(cnames))
    cnames <- paste0("par",1:length(initial))
  
  # Containers
  R            <- runif(nbatch)
  ans          <- matrix(ncol = length(initial), nrow = nbatch,
                         dimnames = list(1:nbatch, cnames))

  # Wrapping function
  f <- function(z) suppressWarnings(fun(z))
  
  # Checking boundaries
  if (length(ub) > 1 && (length(initial) != length(ub)))
    stop("Incorrect length of -ub-")
  
  if (length(lb) > 1 && (length(initial) != length(lb)))
    stop("Incorrect length of -lb-")
  
  if (length(ub) == 1 && (length(initial) != length(ub)))
    ub <- rep(ub, length(initial))

  if (length(lb) == 1 && (length(initial) != length(lb)))
    lb <- rep(lb, length(initial))
  
  if (any(ub <= lb))
    stop("-ub- cannot be <= than -lb-.")
  
  # Checkihg burnins
  if (burnin >= nbatch)
    stop("-burnin- (",burnin,") cannot be >= than -nbatch- (",nbatch,").")

  # Checking thin
  if (thin > nbatch)
    stop("-thin- (",thin,") cannot be > than -nbatch- (",nbatch,").")
  
  theta0 <- initial
  f0     <- f(theta0)
  nskip  <- 0L
  for (i in 1:nbatch) {
    # Step 1. Propose
    theta1 <- normal_prop(theta0, ub, lb, scale)
    f1     <- f(theta1)
    
    # Checking f(theta1) (it must be a number, can be Inf)
    if (is.nan(f1) | is.na(f1)) {
      warning("The output of fun(par) is either -NaN- or -NA-. ",
              "We will skip this iteration for now. ",
              "Check either -fun- or the -lb- and -ub- parameters.")
      
      ans   <- ans[-i,,drop=FALSE]
      nskip <- nskip + 1L
      
      next
    }

    # Minus infinite is interpreted as highly implausible (since proposal was done)
    # within boundaries; hence, we replace f1 with a value -.Machine$double.xmax*inf_eps
    if (f1 == -Inf) {
      f1 <- -.Machine$double.xmax*inf.eps
    } else if (f1 == Inf) {
      f1 <- .Machine$double.xmax*inf.eps
    }
    
    
    # Step 2. Hastings ratio
    r <-
      tryCatch(
        min(1, exp(f1 - f0)),
        error = function(e)
          list(t0 = theta0, t1 = theta1, err = e)
      )
    
    # Checking if there's an error: This should be an error
    if (length(r) > 1 | is.nan(r)) {
      stop("Huston, we have a problem. The value of the Metropolist-Hastings ",
           "ratio is either -NaN- or it returned error.\n", r)
      
      next
    }
    
    # Updating the value
    if (R[i] < r) {
      theta0 <- theta1
      f0     <- f(theta0)
    }
    
    # Storing
    ans[i - nskip,] <- theta0
    
  }
  
  # Thinning the data
  ans <- ans[-c(1:burnin),]
  ans <- ans[(1:nrow(ans) %% thin) == 0, , drop=FALSE]
  
  # Returning an mcmc object from the coda package
  # if the coda package hasn't been loaded, then return a warning
  if (!("package:coda" %in% search()))
    warning("The -coda- package has not been loaded.")
  
  return(
    coda::mcmc(
      ans,
      start = as.integer(rownames(ans)[1]),
      end   = as.integer(rownames(ans)[nrow(ans)]),
      thin = thin
      )
    )
}



#' @rdname mle
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
  if (!length(control$scale))  control$scale  <- .0005
  if (!length(control$ub))     control$ub     <- rep(1 - 1e-20, 5)
  if (!length(control$lb))     control$lb     <- rep(1e-20, 5)
  
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
      
      res <- LogLike(
        Z          = dat$experiment,
        offspring  = dat$offspring,
        noffspring = dat$noffspring,
        psi        = params[1:2] ,
        mu         = params[3:4] ,
        Pi         = c(1 - params[5], params[5]),
        verb_ans   = FALSE
      )$ll + sum(log(priors(params)))
      
      if (is.nan(res) | is.infinite(res)) return(-Inf)
      res
      
    }
  } else {
    function(params) {
      # These are probabilitis, so they should be in the unit sphere.
      # So the probability of these is zero
      if (any(params > 1 | params < 0))
          return(-Inf)
      
      # Checking whether params are fixed or not
      params <- ifelse(fix.params, par0, params)
      
      res <- LogLike(
        Z          = dat$experiment,
        offspring  = dat$offspring,
        noffspring = dat$noffspring,
        psi        = params[1:2] ,
        mu         = params[3:4] ,
        Pi         = c(1 - params[5], params[5]),
        verb_ans   = FALSE
      )$ll
      
      if (is.nan(res) | is.infinite(res)) return(-Inf)
      res
    }
  }
  
  # Naming the parameter estimates
  names(params) <- c("psi0", "psi1", "mu0", "mu1", "Pi")
  
  # Running the MCMC
  ans <- do.call(MCMC, c(list(fun = fun, initial = params), control))
  
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
      par = ans[nrow(ans),],
      hist = ans,
      value = fun(ans[nrow(ans),]),
      ll    = fun(ans[nrow(ans),]),
      fun = fun,
      priors = priors,
      dat = dat,
      par0 = par0,
      fix.params = fix.params,
      method = "mcmc"
    ),
    class = "phylo_mle"
  )
}


