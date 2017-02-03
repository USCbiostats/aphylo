#' Markov Chain Monte Carlo
#' @param fun A function. Returns the log-likelihood
#' @param initial A numeric vector. Initial values of the parameters.
#' @param nbatch Integer scalar. Number of MCMC runs.
#' @param scale Numeric scalar. Step size, \eqn{z\times scale}{z*scale}, where
#' \eqn{z\sim N(0,1)}{z~N(0,1)}.
#' @param burnin Integer scalar. Number of burn-in samples.
#' 
#' @details This function implements MCMC using the Hastings ratio with
#' scaled standard normal propositions for each parameter. For each parameter
#' the transition function is 
#' 
#' \deqn{
#' \theta' = \theta + scale*z
#' }
#' 
#' Where \eqn{z} has standard normal distribution.
#' 
#' @return 
#' \item{batch}{Numeric matrix of size \code{length(initial) x nbatch}.}
#' \item{final}{Numeric vector of length \code{length(initial)}. Final state of the parameters.}
#' \item{initial}{Numeric vector of length \code{length(initial)}. 
#' Initial state of the parameters.}
#' \item{fun}{Objective function.}
#' 
#' @export
#' 
mcmc <- function(
  fun,
  initial, 
  nbatch = 2e4L,
  scale = 1,
  burnin = 1e3L,
  ub = rep(1 - 1e-20, length(initial)),
  lb = rep(1e-20, length(initial))
  ) {
  
  # Containers
  R   <- runif(nbatch)
  ans <- matrix(ncol = length(initial), nrow = nbatch)
  
  # Wrapping function
  f <- function(z) suppressWarnings(fun(z))
  
  # Checking boundaries

  
  theta0 <- initial
  f0     <- f(theta0)
  for (i in 1:nbatch) {
    # Step 1. Propose
    theta1 <- normal_prop(theta0, ub, lb, scale)
    f1     <- f(theta1)
    
    # Check proposition
    if (is.nan(f1) | is.infinite(f1))
      next
    
    # Step 2. Hastings ratio
    r <-
      tryCatch(
        exp(f1 - f0),
        error = function(e)
          list(t0 = theta0, t1 = theta1, err = e)
      )
    
    # Checking if there's an error
    if (length(r) > 1 | is.nan(r) | is.infinite(r)) {
      cat(sprintf("theta0: %f theta1: %f r: %f\n", theta0, theta1, r))
      print(ans[1:i,])
      stop("Ups! Huston, we have a problem.")
      
    }
    
    # Updating the value
    if (R[i] < r) {
      theta0 <- theta1
      f0     <- f(theta0)
    }
    
    # Storing
    ans[i,] <- theta0
    
  }
  
  return(list(
    batch = ans[-c(1:burnin), ],
    final = theta0,
    initial = initial,
    fun   = fun)
    )
  
}


# # Parameters
# set.seed(1231)
# n <- 1e3
# pars <- c(mean = 2.6, sd = 3)
# 
# # Generating data and writing the log likelihood function
# D <- rnorm(n, pars[1], pars[2])
# fun <- function(x) {
#   x <- log(dnorm(D, x[1], x[2]))
#   if (any(is.infinite(x)))
#     return(-Inf)
#   sum(x)
# }
# ans <- mcmc(fun, rep(1,2), nbatch = 2e3)

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
      
      res <- LogLike(dat$experiment, dat$offspring, dat$noffspring, psi, mu, Pi, FALSE)$ll +
        sum(log(priors(params)))
      
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
      
      psi <- params[1:2] 
      mu  <- params[3:4] 
      Pi  <- params[5]
      Pi  <- c(1 - Pi, Pi)
      
      res <- LogLike(dat$experiment, dat$offspring, dat$noffspring, psi, mu, Pi, FALSE)$ll
      
      if (is.nan(res) | is.infinite(res)) return(-Inf)
      res
    }
  }
  
  # ans <- do.call(mcmc::metrop, c(list(obj = fun, initial = params), control))
  ans <- do.call(mcmc, c(list(fun = fun, initial = params), control))
  
  # Working on answer
  env <- new.env()
  environment(fun)    <- env
  environment(dat)    <- env
  environment(par0)   <- env
  environment(fix.params) <- env
  if (length(priors)) 
    environment(priors) <- env
  
  # Naming the parameter estimates
  names(ans$final) <- c("psi0", "psi1", "mu0", "mu1", "Pi")
  
  # Returning
  structure(
    list(
      par = ans$final,
      hist = ans$batch,
      value = fun(ans$final),
      ll    = fun(ans$final),
      fun = fun,
      priors = priors,
      dat = dat,
      par0 = par0,
      fix.params = fix.params,
      method = "mcmc",
      time = ans$time,
      accept = ans$accept
    ),
    class = "phylo_mle"
  )
}