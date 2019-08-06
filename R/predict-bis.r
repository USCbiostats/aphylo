#' Posterior probabilities based on parameter estimates
#' 
#' The function `predict_pre_order` uses a pre-order algorithm to compute the
#' posterior probabilities, whereas the `predict_brute_force` computes posterior
#' probabilities generating all possible cases.
#' 
#' @param atree,x,object Either a tree of class [aphylo] or an object of class [aphylo_estimates]
#' @param params A numeric vector with the corresponding parameters.
#' @param newdata (optional) An aphylo object.
#' 
#' @template parameters
#' @templateVar .psi 1
#' @templateVar .mu 1
#' @templateVar .Pi 1
#' @templateVar .eta 1
#' @param ... Ignored.
#' 
#' @details 
#' The function `predict_brute_force` is only intended for testing. For predictions
#' after estimating the model, see [predict.aphylo_estimates].
#' 
#' @name posterior-probabilities
NULL

#' @rdname posterior-probabilities
#' @return In the case of the `predict` method, a `P` column numeric matrix
#' with values between \eqn{[0,1]} (probabilities).
#' @export
predict.aphylo_estimates <- function(object, newdata = NULL,...) {
  
  # Running prediction function
  pred <- predict_pre_order.aphylo_estimates(object, newdata = newdata, ...)
  
  # No need to check for the class since we are already doing that in
  # predict_preorder.aphylo_estimates
  if (!is.null(newdata)) 
    object$dat <- newdata
  
  # Adding names
  if (is.aphylo(object$dat)) {
    
    dimnames(pred) <- list(
      with(object$dat$tree, c(tip.label, node.label)),
      colnames(object$dat$tip.annotation))
    
  } else {
    
    for (i in seq_along(pred)) {
      dimnames(pred[[i]]) <- list(
        with(object$dat[[i]]$tree, c(tip.label, node.label)),
        colnames(object$dat[[i]]$tip.annotation))
    }
    
  }
  
  
  pred
}

#' @export
#' @rdname posterior-probabilities
predict_pre_order <- function(x, ...) UseMethod("predict_pre_order")

#' @export
#' @rdname posterior-probabilities
predict_pre_order.aphylo_estimates <- function(
  x,
  params  = x$par,
  newdata = NULL,
  ...
  ) {
  
  dots <- list(...)
  
  if (!is.null(newdata)) {
    if (inherits(newdata, "aphylo") | inherits(newdata, "multiAphylo"))
      x$dat <- newdata
    else 
      stop("`newdata` should be of class `aphylo` or `multiAphylo`.",
           call. = FALSE)
  }
  
  if (Ntrees(x) > 1) {
    
    ans <- vector("list", Ntrees(x))
    x.  <- x
    for (t. in seq_along(ans)) {
      x.$dat <- x$dat[[t.]]
      ans[[t.]] <-predict_pre_order(x., params = x$par, ...)
    }
    
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
  p   <- Nann(x)
  ans <- lapply(1:p, function(i) {
    
    # Updating tree (we only need a single function)
    tmpdat <- x$dat[i]
    
    dots$dat      <- tmpdat
    dots$p        <- params
    dots$verb_ans <- TRUE
    
    # The user may change the priors
    if (!length(dots$priors))
      dots$priors <- x$priors
    
    # Computing loglike
    l <- do.call(x$fun, dots)
    
    # Returning posterior probability
    .posterior_prob(
      Pr_postorder = l$Pr[[1]],
      types        = types,
      mu_d         = mu_d,
      mu_s         = mu_s,
      Pi           = Pi,
      pseq         = x$dat$pseq,
      offspring    = x$dat$offspring
      )$posterior
    
  })
  
  do.call(cbind, ans)
  
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
      tree       = x[i],
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

