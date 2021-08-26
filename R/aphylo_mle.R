#' Model estimation using Maximum Likelihood Estimation
#'
#' The function is a wrapper of [stats::optim()].
#' 
#' @template estimates
#' @param method,control,lower,upper Arguments passed to [stats::optim()]. 
#' @family parameter estimation
#' @export
#' @details 
#' The default starting parameters are described in [APHYLO_PARAM_DEFAULT].
#' @examples 
#' 
#' # Using simulated data ------------------------------------------------------
#' set.seed(19)
#' dat <- raphylo(100)
#' dat <- rdrop_annotations(dat, .4)
#' 
#' # Computing Estimating the parameters 
#' ans  <- aphylo_mle(dat ~ psi + mu_d + eta + Pi)
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
#' ans_dbeta <- aphylo_mle(dat ~ psi + mu_d + eta + Pi, priors = mypriors)
#' ans_dbeta
#' @aliases MLE
aphylo_mle <- function(
  model,
  params,
  method            = "L-BFGS-B",
  priors            = function(p) 1, 
  control           = list(),
  lower             = 1e-5,
  upper             = 1 - 1e-5,
  check_informative = getOption("aphylo_informative", FALSE),
  reduced_pseq      = getOption("aphylo_reduce_pseq", TRUE)
) {
  
  # Getting the call
  cl <- match.call()
  
  # Parsing the formula
  env   <- environment(model)
  model <- aphylo_formula(model, params, priors, env = env)
  
  if (any(model$fixed))
    warnings("Fixing parameters is ignored in MLE estimation.")
  
  # If all are 9s, then, there's nothing to do with it.
  if (sum(Nannotated(model$dat)) == 0L) {
    if (check_informative)
      stop("This tree is empty. With no annotations on the tips, no model can be estimated.", call.=FALSE)
    else
      warning("This tree is empty. With no annotations on the tips, no model can be estimated.", call.=FALSE)
  }
  
  # Reducing the peeling sequence
  # This only happens if the eta parameter is not included
  if (length(model$fixed["eta0"]))
    reduced_pseq <- FALSE
  
  # Use the longer the pruning sequence?
  if (!reduced_pseq) {
    old_reduced_pseq       <- model$dat$reduced_pseq
    model$dat$reduced_pseq <- model$dat$pseq
  }
  
  # If the models is uninformative, then it will return with error
  if (check_informative)
    stop_ifuninformative(model$dat$tip.annotation)
  
  # Both available for ABC
  control$fnscale  <- -1
  control$parscale <- rep(1, length(model$params))
  
  # Optimizing
  dat0 <- new_aphylo_pruner(model$dat)
  ans <- do.call(
    stats::optim, 
    c(
      list(
        par      = model$params,
        fn       = model$fun,
        dat      = dat0,
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
    ans$par, model$fun, dat = dat0, priors = priors, verb_ans = FALSE,
    control = control
  )
  
  # Hessian for observed information matrix
  dimnames(hessian) <- list(names(ans$par), names(ans$par))
  
  # Going back to the original reduced pseq
  if (!reduced_pseq) 
    model$dat$reduced_pseq <- old_reduced_pseq
  
  # Returning
  new_aphylo_estimates(
    par         = ans$par,
    hist        = NULL,
    ll          = ans$value,
    counts      = ans$counts,
    convergence = ans$convergence,
    message     = ans$message,
    fun         = model$fun,
    priors      = priors,
    dat         = model$dat,
    par0        = model$params,
    method      = method,
    varcovar    = try_solve(-hessian, tol = 1e-100),
    call        = cl
  )
}
