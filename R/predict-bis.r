#' Posterior probabilities based on parameter estimates
#' 
#' The function `predict_pre_order` uses a pre-order algorithm to compute the
#' posterior probabilities, whereas the `predict_brute_force` computes posterior
#' probabilities generating all possible cases.
#' 
#' @param atree,x,object Either a tree of class [aphylo] or an object of class [aphylo_estimates]
#' @param params A numeric vector with the corresponding parameters.
#' @param newdata (optional) An aphylo object.
#' @param which.tree Integer scalar. Which tree to include in the prediction.
#' 
#' @template loo
#' @template parameters
#' @templateVar .psi 1
#' @templateVar .mu 1
#' @templateVar .Pi 1
#' @templateVar .eta 1
#' @param ... Ignored.
#' 
#' 
#' @details 
#' The function `predict_brute_force` is only intended for testing. For predictions
#' after estimating the model, see [predict.aphylo_estimates].
#' 
#' In the case of the parameter `loo` (leave-one-out), while making tip-level
#' predictions, at each leaf the algorithm will drop annotations regarding that
#' leaf, making its prediction using all the available information except the
#' one include in such leaf.
#' 
#' @name posterior-probabilities
NULL

#' @rdname posterior-probabilities
#' @return In the case of the `predict` method, a `P` column numeric matrix
#' with values between \eqn{[0,1]} (probabilities).
#' @export
predict.aphylo_estimates <- function(
  object,
  which.tree = 1:Ntrees(object),
  ids        = lapply(Ntip(object)[which.tree], seq_len),
  newdata    = NULL,
  loo        = TRUE,
  nsamples   = 1L,
  centiles   = c(.025, .5, .975),
  cl         = NULL,
  ...
  ) {
  
  if (any(Nann(object)[which.tree] > 1L) && (nsamples > 1L))
    stop(
      "Predictions using multiple samples is restricted to a single function ",
      "for now.", call. = FALSE
      )
  
  # Running prediction function
  pred <- predict_pre_order.aphylo_estimates(
    object,
    which.tree = which.tree,
    ids        = ids,
    newdata    = newdata,
    loo        = loo,
    nsamples   = nsamples,
    centiles   = centiles,
    cl         = cl,
    ...
    )
  
  # No need to check for the class since we are already doing that in
  # predict_preorder.aphylo_estimates
  if (!is.null(newdata)) 
    object$dat <- newdata
  
  # Adding names
  if (is.aphylo(object$dat)) {
    
    # Figuring out names
    cnames <- if (nsamples > 1)
      sprintf(
        "%s_%.3f", colnames(object$dat$tip.annotation),
        centiles
        )
    else
      colnames(object$dat$tip.annotation)
    
    dimnames(pred) <- list(
      with(object$dat$tree, c(tip.label, node.label)),
      cnames
      )
    
  } else {
    
    for (i in seq_along(which.tree)) {
      
      # Figuring out names
      cnames <- if (nsamples > 1)
        sprintf(
          "%s_%.3f",
          colnames(object$dat[[which.tree[i]]]$tip.annotation),
          centiles
        )
      else
        colnames(object$dat[[which.tree[i]]]$tip.annotation)
      
      dimnames(pred[[i]]) <- list(
        with(object$dat[[which.tree[i]]]$tree, c(tip.label, node.label)),
        cnames
        )
      
    }
    
  }
  
  pred
}

#' @export
#' @rdname posterior-probabilities
predict_pre_order <- function(x, ...) UseMethod("predict_pre_order")

#' @export
#' @rdname posterior-probabilities
#' @param nsamples Integer scalar. When greater than one, the prediction is done
#' using a random sample from the MCMC chain. This only works if the model was
#' fitted using MCMC, of course.
#' @param ncores,cl Passed to [parallel::makeCluster()].
#' @param ids Integer vector. Ids (positions) of the nodes that need to be
#' predicted.
#' @param centiles Used together with `nsamples`, this indicates the centiles
#' to be computed from the distribution of outcomes.
predict_pre_order.aphylo_estimates <- function(
  x,
  params     = x$par,
  which.tree = 1:Ntrees(x),
  ids        = lapply(Ntip(x)[which.tree], seq_len),
  newdata    = NULL,
  loo        = TRUE,
  nsamples   = 1L,
  centiles   = c(.025, .5, .975),
  ncores     = 1L,
  cl         = NULL,
  ...
  ) {
  
  # Checking that the new data makes sense
  if (!is.null(newdata)) {
    if (inherits(newdata, "aphylo") | inherits(newdata, "multiAphylo"))
      x$dat <- newdata
    else 
      stop("`newdata` should be of class `aphylo` or `multiAphylo`.",
           call. = FALSE)
  }
  
  # Multiple trees are simply passed along the way -----------------------------
  if (is.multiAphylo(x$dat)) {
    
    ans <- vector("list", length(which.tree))
    x.  <- x
    
    for (t. in seq_along(which.tree)) {
      x.$dat <- x$dat[[which.tree[t.]]]
      ans[[t.]] <- predict_pre_order(
        x      = x.,
        ids    = ids[t.],
        params = x$par,
        nsamples = nsamples,
        ncores   = ncores,
        centiles = centiles,
        ...
        )
    }
    
    return(ans)
  }
  
  # In the case of multiple samples --------------------------------------------
  if (nsamples > 1L) {
    
    # Only valid for MCMC
    if (x$method != "mcmc")
      stop(
        "Using the nsamples parameters is only valid for MCMC estimation ",
        "method, this is ", x$method, ".", call. = FALSE
      )
    
    # Valid ranges
    centiles <- sort(centiles)
    if (centiles[1L] < 0 | centiles[length(centiles)] > 1)
      stop("Out of range. When specifying centiles, these should be within the ",
           "[0, 1] range.", call. = FALSE
      )
    
    # Warning b/c of number of samples
    if (nsamples > coda::niter(x$hist) * coda::nchain(x$hist))
      warning("Retrieving more samples than iterations are available.",
              call. = FALSE, immediate. = TRUE
      )
    
    # Sampling parameters
    samples <- do.call(rbind, x$hist) 
    samples <- lapply(
      sample.int(n = nrow(samples), size = nsamples, replace = TRUE), 
      function(i) samples[i, ]
    )
    
    if (ncores > 1L & is.null(cl)) {
      on.exit(tryCatch(parallel::stopCluster(cl), error = function(e) e))
      cl <- parallel::makeCluster(ncores)
    } else if (!is.null(cl)) 
      ncores <- length(cl)
    
    # Calling the prediction function
    ans <- if (ncores > 1L) {
      parallel::parLapply(
        cl,
        X = samples, function(params., x., ids, newdata, loo, ...) {
          aphylo::predict_pre_order(
            x       = x.,
            ids     = ids,
            newdata = newdata,
            loo     = loo,
            params  = params.,
            ...
          )
        },
        x.      = x,
        ids     = ids[1L],
        newdata = newdata,
        loo     = loo,
        ...
      )
    } else {
      
      lapply(
        X = samples, function(params., x, ids, newdata, loo, ...) {
          aphylo::predict_pre_order(
            x       = x,
            ids     = ids,
            newdata = newdata,
            loo     = loo,
            params  = params.,
            ...
          )
        },
        x       = x,
        ids     = ids[1L],
        newdata = newdata,
        loo     = loo,
        ...
      )
      
    }
    
    # Combining results and computing desired centiles
    ans0 <- do.call(cbind, ans)
    ans  <- matrix(nrow = nrow(ans0), ncol = length(centiles))
    ans[ids[[1]], ] <- if (ncores > 1L) {
      t(parallel::parApply(cl, X = ans0[ids[[1]], ], MARGIN = 1L, stats::quantile, prob = centiles))
    } else
      t(apply(ans0[ids[[1]], ], MARGIN = 1L, FUN = stats::quantile, prob = centiles))
    
    # Done! For now, this will be a bit of a mess if we have multiple functions
    return(ans)
    
  }
  
  # Checking parameters
  types <- with(x$dat, c(tip.type, node.type))
  mu_d <- params[c("mu_d0", "mu_d1")]
  if (!("mu_s0" %in% names(params)))
    mu_s <- mu_d
  else
    mu_s <- params[c("mu_s0", "mu_s1")]
  
  if (!("Pi" %in% names(params))) {
    p0 <- mean(types == 0L)
    Pi <-
      p0 * mu_d[1]/(mu_d[1] + mu_d[2]) +
      (1-p0) * mu_s[1]/(mu_s[1] + mu_s[2])
  } else 
    Pi <- params["Pi"]
  
  # Looping through the variables
  dots <- list(...)
  p   <- Nann(x)
  ans <- matrix(nrow = ape::Nnode(x, internal.only = FALSE), ncol = p)
  for (j in 1L:p) {

    # Updating tree (we only need a single function)
    tmpdat <- x$dat[,j]
    
    dots$dat      <- tmpdat
    dots$p        <- params
    dots$verb_ans <- TRUE
    
    # The user may change the priors
    if (!length(dots$priors))
      dots$priors <- x$priors
    
    # Computing loglike
    if (!loo) {
    
      l <- do.call(x$fun, dots)
      
      # Returning posterior probability
      ans[, j] <- .posterior_prob(
        Pr_postorder = l$Pr[[1]],
        types        = types,
        mu_d         = mu_d,
        mu_s         = mu_s,
        Pi           = Pi,
        pseq         = x$dat$pseq,
        offspring    = x$dat$offspring
        )$posterior
      
      next
      
    }
    
    dots$dat <- new_aphylo_pruner(tmpdat)
    for (i in intersect(1L:Ntip(x), ids[[1L]])) {
      
      # Setting that annotation to Missing (9)
      Tree_set_ann(dots$dat, i - 1L, 0L, 9L)
      l <- do.call(x$fun, dots)
      
      ans[i, j] <- .posterior_prob(
        Pr_postorder = l$Pr[[1L]],
        types        = types,
        mu_d         = mu_d,
        mu_s         = mu_s,
        Pi           = Pi,
        pseq         = x$dat$pseq,
        offspring    = x$dat$offspring
      )$posterior[i,]
      
      # Original value
      Tree_set_ann(dots$dat, i - 1L, 0L, x$dat$tip.annotation[i, j])
      
    }
    
    # Filling the rest of the tree
    l <- do.call(x$fun, dots)
    last_set <- intersect((i + 1L):nrow(ans), ids[[1L]])
    
    if (length(last_set))
      ans[last_set, j] <- .posterior_prob(
        Pr_postorder = l$Pr[[1L]],
        types        = types,
        mu_d         = mu_d,
        mu_s         = mu_s,
        Pi           = Pi,
        pseq         = x$dat$pseq,
        offspring    = x$dat$offspring
      )$posterior[last_set, ]
    
  }
  
  ans
  
}

#' @rdname posterior-probabilities
#' @export
predict_pre_order.aphylo <- function(x, psi, mu_d, mu_s, eta, Pi, ...) {
  
  if (Ntrees(x) > 1) {
    
    ans <- vector("list", Ntrees(x))
    # x.  <- x
    for (t. in seq_along(ans)) {
      x.$dat <- x$dat[[t.]]
      ans[[t.]] <- predict_pre_order(x., psi, mu_d, mu_s, eta, Pi, ...)
    }
    
    return(ans)
      
  }
  
  p <- ncol(x$tip.annotation)
  ans <- lapply(1:p, function(i) {
    
    # Computing loglike
    l <- LogLike(
      tree       = x[, i],
      psi        = psi,
      mu_d       = mu_d,
      mu_s       = mu_s,
      eta        = eta,
      Pi         = Pi,
      verb_ans   = TRUE,
      check_dims = FALSE
    )
    
    # Returning posterior probability
    .posterior_prob(
      Pr_postorder = l$Pr[[1L]],
      types        = with(x, c(tip.type, node.type)),
      mu_d         = mu_d,
      mu_s         = mu_s,
      Pi           = Pi,
      pseq         = x$pseq,
      offspring    = x$offspring
      )$posterior
    
  })
    
  do.call(cbind, ans)
  
}

#' @rdname posterior-probabilities
#' @export
predict_brute_force <- function(atree, psi, mu_d, mu_s, Pi) {
  
  # Should be aphylo
  if (!inherits(atree, "aphylo"))
    stop("`atree` must be of class `aphylo` (it is of class ", class(atree), ".")
  
  if (length(atree$offspring) > 7)
    stop("In the case of the brute-force calculations, trees with more than 7 nodes becomes burdensome.")
  
  # Coercing into the true class
  tree <- as.phylo(atree)
  Pi   <- c(1 - Pi, Pi)
  
  # Generating a states matrix
  states <- do.call(expand.grid, rep(list(c(0,1)), ape::Nnode(tree) + ape::Ntip(tree)))
  colnames(states) <- paste0("node", 1:(ncol(states)))
  
  # Computing matrix of probabilities --------------------------------------------
  # 2^(ntips + nnodes) 
  PSI <- prob_mat(psi)
  MU  <- list(prob_mat(mu_d), prob_mat(mu_s))
  
  # For
  Pr <- Pi[states[, ape::Ntip(tree) + 1] + 1]
  types <- with(atree, c(tip.type, node.type))
  for (i in 1:nrow(tree$edge)) {
    e <- tree$edge[i,]
    Pr <- Pr *
      MU[[ types[e[1]] + 1L ]][cbind(states[, e[1]], states[, e[2]]) + 1]
  }
  
  # Computing for each possible annotation of the tips
  # 2^Ntip
  Pr <- matrix(rep(Pr, 3), nrow=length(Pr), ncol = 2^ape::Ntip(tree))
  states_tip <- do.call(expand.grid, rep(list(c(0,1)), ape::Ntip(tree)))
  
  for (i in 1:nrow(states_tip)) {
    
    for (j in 1:ape::Ntip(tree)) {
      Pr[, i] <- Pr[, i] *
        PSI[cbind(states[, j] + 1, states_tip[i, j] + 1)]
    }
    
  }
  
  # Computing posterior probabilities
  posterior <- vector("numeric", ncol(states))
  for (i in 1:ncol(states)) {
    state_col <- which(apply(states_tip, 1, function(x) all(x == atree$tip.annotation)))
    state_rows <- which(states[,i] == 1)
    
    posterior[i] <- sum(Pr[state_rows, state_col])/sum(Pr[, state_col])
  }
  
  
  # Results
  list(
    Pr  = Pr,
    row = states,
    col = states_tip,
    posterior = posterior
  )
  
}

