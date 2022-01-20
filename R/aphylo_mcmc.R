#' @export
#' @rdname aphylo_mcmc
#' @details `APHYLO_DEFAULT_MCMC_CONTROL` lists the default values for the MCMC
#' estimation:
#' - `nsteps`: `1e4L`
#' - `burnin`: `5e3L`
#' - `thin` : `10L`
#' - `nchains` : `2L`
#' - `multicore` : `FALSE`
#' - `conv_checker` : `fmcmc::convergence_auto(5e3)`
#' 
#' For more information about the MCMC estimation process, see [fmcmc::MCMC()].
APHYLO_DEFAULT_MCMC_CONTROL <- list(
  nsteps    = 1e4L,
  burnin    = 5e3L,
  thin      = 10L,
  nchains   = 2L,
  multicore = FALSE,
  conv_checker = fmcmc::convergence_auto(5e3)
)

#' Model estimation using Markov Chain Monte Carlo
#'
#' The function is a wrapper of [fmcmc::MCMC()].
#'
#' @template estimates
#' @param control A list with parameters for the optimization method (see
#' details).
#' @details 
#' Methods [base::print()], [base::summary()], [stats::coef], [stats::window()],
#' [stats::vcov()], [stats::logLik()], [predict()][predict.aphylo_estimates()],
#' and the various ways to query features of the trees via [Ntip()][ape::Ntip()]
#' are available post estimation.
#' 
#' @family parameter estimation
#' @export
#' @examples 
#' # Using the MCMC ------------------------------------------------------------
#' 
#' \donttest{
#' 
#' set.seed(1233)
#' # Simulating a tree
#' tree <- sim_tree(200)
#' 
#' # Simulating functions
#' atree <- raphylo(
#'   tree = tree,
#'   psi  = c(.01, .03),
#'   mu_d = c(.05, .02),
#'   Pi   = .5
#' )
#' 
#' # Running the MCMC
#' set.seed(1231)
#' 
#' ans_mcmc <- aphylo_mcmc(
#'   dat ~ mu_d + psi + eta + Pi,
#'   control = list(nsteps = 2e5, burnin=1000, thin=200)
#' )
#' }
#' @aliases MCMC
aphylo_mcmc <- function(
  model,
  params,
  priors            = uprior(),
  control           = list(),
  check_informative = getOption("aphylo_informative", FALSE),
  reduced_pseq      = getOption("aphylo_reduce_pseq", TRUE)
) {
  
  # Getting the call
  cl <- match.call()
  
  # Parsing the formula
  env   <- environment(model)
  model <- aphylo_formula(model, params, priors, env = env)
  
  # If all are 9s, then, there's nothing to do with it.
  if (sum(Nannotated(model$dat)) == 0L) {
    if (check_informative)
      stop("This tree is empty. With no annotations on the tips, no model can be estimated.", call.=FALSE)
    else
      warning("This tree is empty. With no annotations on the tips, no model can be estimated.", call.=FALSE)
  }
  
  # Reducing the peeling sequence
  # This only happens if the eta parameter is not included
  if ("eta0" %in% names(model$fixed))
    reduced_pseq <- FALSE
  
  # Use the longer the pruning sequence?
  if (!reduced_pseq) {
    old_reduced_pseq       <- model$dat$reduced_pseq
    model$dat$reduced_pseq <- model$dat$pseq
  }
  
  # Checking control
  for (n in names(APHYLO_DEFAULT_MCMC_CONTROL)) {
    if (!length(control[[n]]) && !(n %in% names(control)))
      control[[n]] <- APHYLO_DEFAULT_MCMC_CONTROL[[n]]
  }
  
  if (!("kernel" %in% names(control)))
    control$kernel <- fmcmc::kernel_am(
      ub   = 1,
      lb   = 0,
      freq = 1L
    )
  
  # If the models is uninformative, then it will return with error
  if (check_informative)
    stop_ifuninformative(model$dat$tip.annotation)
  
  if ("multicore" %in% names(control) && control$multicore) {
    
    # Initalizing the cluster
    cl_object <- parallel::makePSOCKcluster(control$nchains)
    on.exit(parallel::stopCluster(cl_object))
    
    # Setting up the package
    parallel::clusterEvalQ(cl_object, library(aphylo))
    parallel::clusterExport(cl_object, c("model", "priors"), envir = environment())
    parallel::clusterEvalQ(cl_object, {
      dat00 <- aphylo::new_aphylo_pruner(model$dat)
    })
    parallel::clusterEvalQ(cl_object, {
      f_pll <- model$fun
      model$fun <- function(p, dat, priors, verb_ans, ...) {
        f_pll(p, dat = dat00, priors = priors, verb_ans = verb_ans)
      }
    })
    
    # Appending to the set of controls
    control$cl <- cl_object
    dat0 <- model$dat
    
  } else 
    dat0 <- new_aphylo_pruner(model$dat)
  
  # Running the MCMC
  ans <- do.call(
    fmcmc::MCMC, 
    c(
      list(
        fun      = model$fun,
        dat      = dat0,
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
  
  ll <- model$fun(
    p        = par,
    dat      = dat0,
    priors   = priors,
    verb_ans = FALSE
  )
  
  # Use the longer the pruning sequence?
  if (!reduced_pseq) 
    model$dat$reduced_pseq <- old_reduced_pseq
  
  # Returning
  new_aphylo_estimates(
    par         = par,
    hist        = ans,
    ll          = ll,
    counts      = as.integer(rownames(ans[[1]])[coda::niter(ans)]),
    convergence = NA,
    message     = NA,
    fun         = model$fun,
    priors      = priors,
    dat         = model$dat,
    par0        = model$params,
    method      = "mcmc",
    varcovar    = stats::var(do.call(rbind, ans)),
    call        = cl
  )
}

# @rdname aphylo_mcmc
#' @export
window.aphylo_estimates <- function(x, ...) {
  
  if (x$method != "mcmc") {
    warning("No window method for aphylo_estimates using MLE.", call. = FALSE)
    return(x)
  }
  
  x$hist <- window(x$hist, ...)
  x$par  <- colMeans(do.call(rbind, x$hist))
  x$ll   <- x$fun(
    p = x$par, dat = x$dat, priors = x$priors, verb_ans = FALSE
  )
  x$varcovar <- stats::var(do.call(rbind, x$hist))
  x
}
