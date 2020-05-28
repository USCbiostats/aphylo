library(aphylo)

#' Fit a Hierarchical aphylo model
#' @param model A model in the form of [aphylo_formula].
#' @param params Set of initial parameters, passed to [aphylo_formula].
#' @param ... Further arguments passed to [fmcmc::MCMC()].
#' @param params0 Starting parameters.
#' @param env Environment where to evaluate `model`.
#' @noRd
aphylo_hier <- function(
  model,
  params,
  ...,
  params0      = NULL,
  hyper_params = NULL,
  env          = parent.frame()
  ) {
  
  # Building the likelihoods
  LHS <- eval(model[[2]], envir = env)
  
  # Checks
  N <- Ntrees(LHS)
  if (N < 2L)
    stop(
      "Hierarchical model should have at least two trees.",
      call. = FALSE
      )
  
  functions <- lapply(LHS, function(d.) {
    model. <- stats::update.formula(model, d. ~ .)
    aphylo_formula(model., env = environment(), params = params)
  })
  
  LHS <- lapply(functions, function(f.) new_aphylo_pruner(f.$dat))
  
  Npar        <- length(functions[[1]]$params)
  par_names0  <- names(functions[[1]]$params)
  par_names   <- lapply(1:N, function(i) sprintf("%s_tree%03i", par_names0, i))
  functions   <- lapply(functions, "[[", "fun")
  
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
    
    sum(likelihoods) 
    
  }
  
  # Figuring out the right prior
  alpha_names <- sprintf("alpha_%s", par_names0)
  beta_names  <- sprintf("beta_%s", par_names0)
  hprior <- function(x, alpha, beta) {
    sum(dbeta(x, shape1 = alpha, shape2 = beta, log = TRUE))
  }
  
  # Initial parameters
  # debug(joint)
  if (is.null(params0))
    params0 <- structure(
      c(rep(params, N), rep(1, Npar * 2)),
      names = c(sprintf("%s_tree%03i", rep(par_names0, N), 1:N), alpha_names, beta_names)
    )
  
  if (!is.null(hyper_params)) {
    params0[alpha_names] <- hyper_params[par_names0, "alpha"]
    params0[beta_names] <- hyper_params[par_names0, "beta"]
  }
    
  
  
  # Adding random noise
  params0   <- rbind(params0, params0)
  params0[] <- jitter(params0, 20)
  params0[params0 < 0] <- 1e-5
  params0[params0 > 1] <- 1 - 1e-5
    
  fmcmc::MCMC(
    initial = params0,
    fun     = joint,
    data.   = LHS,
    hprior  = hprior, 
    ...
    )
  
  
}

# Generating the data ----------------------------------------------------------
set.seed(123)

mu_d <- c(.9, .5)
mu_s <- c(.1, .1)
psi  <- c(.1, .2)
Pi   <- .2

trees <- replicate(
  20, {
    rdrop_annotations(
      raphylo(50, psi = psi, mu_d=mu_d, mu_s = mu_s, Pi = Pi),
      .5
      )
    },
  simplify = FALSE
  )
trees <- do.call(c, trees)

model <- lapply(trees, function(t.) aphylo_formula(t. ~ psi + mu_d + mu_s + Pi))

# Fitting the model ------------------------------------------------------------

# Getting a reference
ans_mle <- aphylo_mle(trees ~ psi + mu_d + mu_s + Pi)

# Constraint matrix
constr <- kronecker(diag(Ntrees(trees) + 2), matrix(TRUE, nrow = 7, ncol = 7))
constr[nrow(constr) - 1:14 + 1,] <- TRUE
constr <- lower.tri(constr)*constr
diag(constr) <- TRUE
Matrix::image(as(constr, "dgCMatrix"))
Matrix::image(as(crossprod(t(constr)), "dgCMatrix"))


k_ram <- fmcmc::kernel_ram(
  lb     = .00000001,
  warmup = 1000,
  eps    = .0001,
  ub     = c(rep(.9999999, 7 * Ntrees(trees)), rep(100, 7 * 2)),
  fixed  = c(rep(FALSE, 7*Ntrees(trees)), rep(TRUE, 7*2)) #,
  # constr = constr
)

k_am <- fmcmc::kernel_am(
  lb = .00000001,
  warmup = 1000,
  eps    = .01,
  fixed  = c(rep(FALSE, 7*Ntrees(trees)), rep(TRUE, 7*2)),
  ub     = c(rep(.9999999, 7 * Ntrees(trees)), rep(100, 7 * 2)), freq = 50
)

# Estimating empirical bayes using map
ALPHAS <- c(psi0 = 2, psi1 = 2, mu_d0 = 18, mu_d1 = 10, mu_s0 = 2, mu_s1 = 2, Pi = 2)
BETAS  <- c(psi0 = 18, psi1 = 18, mu_d0 = 2, mu_d1 = 10, mu_s0 = 18, mu_s1 = 18, Pi = 18)

bpriors <- bprior(shape1 = ALPHAS, shape2 = BETAS)

map_estimates <- vector("list", Ntrees(trees))
for (i in seq_along(trees)) {
  map_estimates[[i]] <- aphylo_mle(trees[[i]] ~ psi + mu_d + mu_s + Pi, priors = bpriors)
  message("Tree ",i, " of ", Ntrees(trees), " done.")
}

map_estimates <- do.call(rbind, lapply(map_estimates, coef))
ab <- vector("list", ncol(map_estimates))
for (i in 1:ncol(map_estimates)) {
  ab[[i]] <- MASS::fitdistr(
    map_estimates[, i], densfun = "beta",
    start = list(shape1 = ALPHAS[i], shape2 = BETAS[i]),
    lower = .001
    )
}


op <- par(mfrow = c(3, 3))
for (i in 1:length(ab))  {
  curve(dbeta(x, ab[[i]]$estimate[1], ab[[i]]$estimate[2]),
        main = names(ALPHAS)[i])
  
  curve(dbeta(x, ALPHAS[i], BETAS[i]),
        main = names(ALPHAS)[i], add = TRUE, col = "red")
}

par(op)

# Compiling priors so we can use them to fix the hyperpriors
alpha_and_beta <- structure(
  t(sapply(ab, coef)),
  dimnames = list(names(ALPHAS), c("alpha", "beta"))
)



set.seed(123)
# debug(k_ram$proposal)
ans0 <- aphylo_hier(
  trees ~ psi + mu_d + mu_s + Pi,
  params       = structure(runif(7), names = names(coef(ans_mle))),
  nsteps       = 1e5,
  kernel       = k_ram,
  thin         = 1,
  nchains      = 2,
  conv_checker = fmcmc::convergence_gelman(freq = 2e3),
  hyper_params = alpha_and_beta
  )

saveRDS(ans0, file = "playground/hierarchical.rds")

# Finding parameters using MLE
if (FALSE) {
  f <- fmcmc::last_("fun")
  ans_mle2 <- optim(
    par    = tail(ans0, 0)[[1]][1,],
    fn     = f,
    data.  = fmcmc::last_("data."),
    hprior = fmcmc::last_("hprior"),
    lower  = k_ram$lb,
    upper  = k_ram$ub,
    control = list(fnscale = -1),
    method = "L-BFGS-B",
    hessian = TRUE
    )
  
  fmcmc::last_("fun")(
    par   = colMeans(do.call(rbind,window(ans0, start = 1000))),
    data. = fmcmc::last_("data."), 
    hprior = fmcmc::last_("hprior")
  )
  
  fmcmc::last_("fun")(
    par    = ans_mle2$par,
    data.  = fmcmc::last_("data."), 
    hprior = fmcmc::last_("hprior")
  )
}

plotter_1 <- function(x, burnin = NULL, ...) {
  
  if (is.null(burnin))
    burnin <- floor((end(x) - start(x))/2)
  
  cnames <- if (is.list(x)) colnames(x[[1]])
  else  colnames(x)
  
  x <- window(x, start = burnin)
  op <- par(mfcol = c(4, 2))
  for (i in 1:4) {
    coda::traceplot(x[,grepl("mu_d0", cnames)][,1:4 + (i-1)*4], ...)#, ylim = c(0,1))
    coda::traceplot(x[,grepl("mu_s0", cnames)][,1:4 + (i-1)*4], ...)#, ylim = c(0,1))
  }
  par(op)
}

plotter_1(ans0, burnin = 1, ylim = c(0,1))

# Plot of distribution of averages

plotter_2 <- function(x, burnin = NULL) {
  
  if (is.null(burnin))
    burnin <- floor((end(x) - start(x))/2)
  
  vars <- c("mu_d0", "mu_d1", "mu_s0", "mu_s1", "Pi")
  op <- par(mfrow = c(3,2))
  for (v in vars) {
    
    # Getting the corresponding alpha
    a_names <- window(x, start = burnin)[,paste0("alpha_", v)]
    b_names <- window(x, start = burnin)[,paste0("beta_", v)]
    
    if (is.list(a_names)) {
      for (i in seq_along(a_names))
        a_names[[i]] <- a_names[[i]]/(a_names[[i]] + b_names[[i]])
    } else
      a_names <- a_names/(a_names + b_names)
    
    
    # plot(density(a_names/(a_names + b_names)),
    #      # xlim = c(0,1),
    #      main = v)
    
    coda::traceplot(a_names, main = v, type="l")
    # abline(h=quantile(a_names, .5), lwd=2, lty=2, col="red")
    
  }
  par(op)
}

plotter_2(ans0, burnin = 1)

stop("Up to here, bud")

set.seed(12315)
k_ram <- fmcmc::kernel_ram(
  lb     = .00000001,
  warmup = 1000,
  eps    = .000001,
  ub     = c(rep(.9999999, 7 * Ntrees(trees)), rep(100, 7 * 2)) #,
  # constr = constr
)

ans1 <- aphylo_hier(
  trees ~ psi + mu_d + mu_s + Pi,
  params = coef(ans_mle),
  nsteps = 5e3,
  kernel = k_ram,
  params0 = rbind(ans_mle2$par)
)

plotter_1(ans1)

# # Plot of distribution of averages
# op <- par(mfrow = c(3,2))
# for (v in vars) {
#   
#   # Getting the corresponding alpha
#   a_names <- window(ans0, start = 1000)[, paste0(v, )]
#   b_names <- window(ans0, start = 1000)[, paste0(v, )]
#   
#   plot(density(a_names/(a_names + b_names)),
#        xlim = c(0,1),
#        main = v)
#   abline(v=quantile(a_names/(a_names + b_names), .5), lwd=2, lty=2, col="red")
#   
# }
# par(op)

# plot(ans0[,c("alpha", "beta"),drop=FALSE])