#' Fit a Hierarchical aphylo model
#' @template estimates
#' @param params0 Starting parameters.
#' @param env Environment where to evaluate `model`.
#' @param ...,multicore,nchains passed to [fmcmc::MCMC]
#' @param classes An integer vector of length equal to the number of trees in
#' the model.
#' @param hyper_params If not specified, the function sets as initial hyper
#' parameters equal to 10.
#' @param verbose Logical scalar. When `TRUE` prints information.
#' @param use_optim Logical, When true uses [stats::optim] as a starting point.
#' @family parameter estimation
#' @details The parameters `priors`, `check_informative`, and `reduced_pseq`
#' are silently ignored in this function.
#' @export
aphylo_hier <- function(
  model,
  params,
  classes,
  ...,
  params0      = NULL,
  hyper_params = NULL,
  env          = parent.frame(),
  verbose      = TRUE,
  multicore    = FALSE,
  nchains      = 1L,
  use_optim    = TRUE,
  priors       = NULL,
  check_informative = NULL,
  reduced_pseq = NULL
  ) {
  
  # Retrieving the trees
  LHS <- eval(model[[2]], envir = env)
  N <- Ntrees(LHS)
  
  # checks ---------------------------------------------------------------------
  if (N < 2L)
    stop("Hierarchical model should have at least two trees.", call. = FALSE)
  
  if (length(classes) != N)
    stop(
      "The -classes- argument should have the same length as the number of trees.",
      call. = FALSE
      )
  
  Nclasses <- length(unique(classes))
  if (length(unique(classes)) == 1)
    warning("Using a single class", call. = FALSE, immediate. = TRUE)
  
  # Building the likelihoods
  formulae <- lapply(LHS, function(d.) {
    model. <- stats::update.formula(model, d. ~ .)
    aphylo_formula(model., env = environment(), params = params)
  })
  
  data. <- lapply(formulae, function(f.) new_aphylo_pruner(f.$dat))
  
  # Building the parameter names (# pars x # classes + # parameters * 2), and
  # the right likelihood function.
  Npar        <- length(formulae[[1]]$params)
  par_names0  <- names(formulae[[1]]$params)
  par_names   <- lapply(1:N, function(i) sprintf("%s_class%03i", par_names0, classes[i]))
  functions   <- lapply(formulae, "[[", "fun")
  
  likelihoods <- double(N)
  
  joint <- function(par, data., hprior) {
    
    for (i in 1L:N)
      likelihoods[i] <- functions[[i]](
        p        = structure(par[par_names[[i]]], names = par_names0),
        dat      = data.[[i]],
        verb_ans = FALSE,
        priors   = function(p) 1
        ) +
        hprior(par[ par_names[[i]] ], alpha = par[alpha_names], beta = par[beta_names])
    
    ans <- sum(likelihoods) 
    if (is.finite(ans))
      return(ans)
    
    return(-.Machine$double.xmax * 1e-8)
    
  }
  
  # Hyper-prior parameters
  alpha_names <- sprintf("alpha_%s", par_names0)
  beta_names  <- sprintf("beta_%s", par_names0)
  hprior <- function(x, alpha, beta) {
    sum(dbeta(x, shape1 = alpha, shape2 = beta, log = TRUE))
  }
  
  # Initial parameters
  if (is.null(params0))
    params0 <- structure(
      c(rep(params, Nclasses), rep(10, Npar * 2)),
      names = c(sprintf("%s_class%03i", rep(par_names0, Nclasses), 1:Nclasses), alpha_names, beta_names)
    )
  
  if (!is.null(hyper_params)) {
    params0[alpha_names] <- hyper_params[par_names0, "alpha"]
    params0[beta_names] <- hyper_params[par_names0, "beta"]
  }
  
  # # Adding random noise
  # nchains <- list(...)$nchains
  # params0 <- t(replicate(nchains, params0))
  # 
  # are_within1 <- which(params0[1,] < 1)
  # 
  # params0[] <- jitter(params0, 30)
  # 
  # params0[,are_within1][params0[,are_within1] > 1] <- 1 - .1
  # params0[params0 < 0] <- .1
  
  # start_point <- ABCoptim::abc_optim(
  #   par     = params0[1,],
  #   fn      = joint,
  #   hprior  = hprior,
  #   data.   = data.,
  #   fnscale = -1,
  #   lb      = rep(.00001, ncol(params0)),
  #   ub      = k_ram$ub
  #   )
  
  
  if (use_optim) {
    # Finding suitable starting point
    if (verbose)
      message("Trying to maximize using L-BFSG-B...")
    start_point <- stats::optim(
      par     = params0,
      fn      = joint,
      hprior  = hprior,
      data.   = data.,
      method  = "L-BFGS-B",
      control = list(fnscale = -1),
      lower   = .0001,
      upper   = c(rep(.9999, Npar * Nclasses), rep(2e3, Npar * 2))
    )
    
    if (verbose)
      message(sprintf(
        "Optimization complete.\n  convergence: %i\n counts: %i\n message: %s",
        start_point$convergence,
        start_point$counts,
        start_point$message
        ))
  } else
    start_point <- list(par = params0)
  
  if (verbose)
    message("Starting MCMC...")
  
  # Preparing the cluster
  if (multicore) {
    
    cl <- parallel::makePSOCKcluster(nchains)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(
      cl,
      c(
        "LHS", "N", "par_names", "par_names0", "alpha_names", "beta_names",
        "likelihoods", "functions"
        ),
      envir = environment()
      )
    
    parallel::clusterEvalQ(cl, {
      library(aphylo)
      library(fmcmc)
      data. <- lapply(LHS, new_aphylo_pruner)
      rm(LHS) # Not needed any-longer
    })
    
    # Updating the joint
    formals(joint) <- list(par = NULL, hprior = NULL)
    environment(joint) <- .GlobalEnv
    
    ans <- fmcmc::MCMC(
      initial   = start_point$par,
      fun       = joint,
      hprior    = hprior, 
      cl        = cl,
      multicore = TRUE, 
      nchains   = nchains,
      progress  = FALSE,
      ...
    )
    
    formals(joint)       <- list(par = NULL, hprior = NULL, data. = NULL)
    environment(joint)   <- environment()
    assign("joint", joint, envir = fmcmc::LAST_MCMC)
    assign("data.", data., envir = fmcmc::LAST_MCMC)
    
  } else {
    
    ans <- fmcmc::MCMC(
      initial   = start_point$par,
      fun       = joint,
      data.     = data.,
      hprior    = hprior, 
      cl        = NULL,
      multicore = FALSE, 
      nchains   = nchains,
      progress  = verbose,
      ...
    )
    
    formals(joint)       <- list(par = NULL, hprior = NULL, data. = NULL)
    environment(joint)   <- environment()
    assign("joint", joint, envir = fmcmc::LAST_MCMC)
    assign("data.", data., envir = fmcmc::LAST_MCMC)
  }
  
  # We treat all chains as mcmc.list
  if (!inherits(ans, "mcmc.list"))
    ans <- coda::mcmc.list(ans)
  
  sol <- colMeans(do.call(rbind, ans))
  new_aphylo_estimates(
    par         = colMeans(do.call(rbind, ans)),
    hist        = ans,
    ll          = joint(sol, data. = data., hprior = hprior),
    counts      = coda::niter(ans),
    convergence = NA,
    message     = NA,
    fun         = joint,
    priors      = function(i) 1,
    dat         = LHS,
    par0        = params0,
    method      = "mcmc",
    varcovar    = stats::var(do.call(rbind, ans)),
    call        = match.call()
  )
 
  # return(ans)
     
}
