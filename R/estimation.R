
#' Try to compute the inverse of the matrix
#' @return If it fails, then it returns a matrix of size n*n with NAs.
#' @noRd
try_solve <- function(x, ...) {
  
  ans <- tryCatch(MASS::ginv(x, ...), error = function(e) e)
  
  # If it is an error
  if (inherits(ans, "error")) {
    warning("The algorithm did not converge. Cannot find the inverse of the hessian.",
            call. = FALSE)
    return(matrix(ncol = length(x), nrow=length(x)))
  }
    
  
  ans
}

#' Objects of class `aphylo_estimates`
#'
#' The model fitting of annotated phylogenetic trees can be done using either
#' MLE via [aphylo_mle()] or MCMC via [aphylo_mcmc()]. This section describes
#' the object of class `aphylo_estimates` that these functions generate and
#' the post estimation methods/functions that can be used.
#' @param x,object Depending of the method, an object of class `aphylo_estimates`.
#' @param ... Further arguments passed to the corresponding method.
#' @return 
#' Objects of class `aphylo_estimates` are a list withh the following elements:
#' \item{par}{A numeric vector of length 5 with the solution.}
#' \item{hist}{A numeric matrix of size `counts*5` with the solution path
#' (length 2 if used `optim` as the intermediate steps are not available to the
#' user). In the case of `aphylo_mcmc`, `hist` is an object of class
#' [coda::mcmc.list()].}
#' \item{ll}{A numeric scalar with the value of `fun(par, dat)`. The value of the log likelihood.}
#' \item{counts}{Integer scalar number of steps/batch performed.}
#' \item{convergence}{Integer scalar. Equal to 0 if `optim` converged. See `optim`.}
#' \item{message}{Character scalar. See `optim`.}
#' \item{fun}{A function (the objective function).}
#' \item{priors}{If specified, the function `priors` passed to the method.}
#' \item{dat}{The data `dat` provided to the function.}
#' \item{par0}{A numeric vector of length 5 with the initial parameters.}
#' \item{method}{Character scalar with the name of the method used.}
#' \item{varcovar}{A matrix of size 5*5. The estimated covariance matrix.}
#' 
#' @name aphylo_estimates
#' @examples 
#' set.seed(7881)
#' atree <- raphylo(40, P = 2)
#' res   <- aphylo_mcmc(atree ~ mu_d + mu_s + Pi)
#' 
#' print(res)
#' coef(res)
#' vcov(res)
#' plot(res)
NULL

new_aphylo_estimates <- function(
  par,
  hist,
  ll,
  counts,
  convergence,
  message,
  fun,
  priors,
  dat,
  par0,
  method,
  varcovar,
  call
) {
  
  structure(
    list(
      par         = par,
      hist        = hist,
      ll          = ll,
      counts      = counts,
      convergence = convergence,
      message     = message,
      fun         = fun,
      priors      = priors,
      dat         = dat,
      par0        = par0,
      method      = method,
      varcovar    = varcovar,
      call        = call
    ),
    class = "aphylo_estimates"
  )
}

#' Stop if the model is uninformative
#' @noRd
stop_ifuninformative <- function(tip.annotation) {
  tab <- tabulate(tip.annotation + 1, 10)[c(1,2,10)]
  # tab <- fast_table_using_labels(tip.annotation, c(0L, 1L))
  if (tab[1L] == 0L | tab[2] == 0L)
    stop("The model is uninformative (either there's only 0s or 1s).", call. = FALSE)
}

#' @export
#' @rdname aphylo_estimates
print.aphylo_estimates <- function(x, ...) {
  # Function to print a bar with variable width
  catbar <- function() paste0(rep("-",options()$width), collapse="")
  
  sderrors   <- structure(sqrt(diag(x$varcovar)), names = names(coef(x)))
  
  ans <- sprintf(
    "\n # of Leafs: %i\n # of Functions %i\n # of Trees: %i\n",
    sum(Ntip(x)), sum(Nann(x)), Ntrees(x)
    )
  ans <- c(ans, sprintf("\n %-6s  %6s  %6s", "", "Estimate", "Std. Err."))
  for (p in names(x$par)) {
    ans <- c(
      ans,
      with(x, sprintf("\n %-6s  %6.4f    %6.4f", p, par[p], sderrors[p]))
      )
  }

  with(x, {
    cat(
      sep = "",
      "\nESTIMATION OF ANNOTATED PHYLOGENETIC TREE\n",
      "\n Call: ", paste(deparse(x$call), sep="\n", collapse="\n"), 
      sprintf(
        "\n LogLik%s: %-9.4f\n Method used: %s (%i steps)",
        if (prod(priors(par)) != 1) " (unnormalized)" else "",
        ll, method, x$counts),
      if (method == "mcmc") 
        NULL
      else
        sprintf("\n convergence: %i (see ?optim)", convergence)
      ,
      ans, "\n\n"
      )
    
  })
  
  invisible(x)
}

#' @export
#' @rdname aphylo_estimates
coef.aphylo_estimates <- function(object, ...) {
  object$par
}

#' @export
#' @rdname aphylo_estimates
vcov.aphylo_estimates <- function(object, ...) {
  object$varcovar
}



#' @export
logLik.aphylo_estimates <- function(object, ...) {
  
  ans <- with(object, fun(par, priors = priors, dat = dat, verb_ans = TRUE))
  
  structure(
    .Data = ans$ll,
    class = "logLik",
    df    = length(object$par),
    Pr    = ans$Pr
  )
  
}

#' @rdname aphylo_estimates
#' @param which.tree Integer scalar. Which tree to plot. 
#' @details The plot method for the object of class `aphylo_estimates` plots
#' the original tree with the predicted annotations.
#' @param y Ignored.
#' @template loo
#' @param ids,nsamples,ncores,centiles,cl passed to [predict.aphylo_estimates()]
#' @return The plot method for `aphylo_estimates` returns the selected tree
#' (`which.tree`) with predicted annotations, also of class [aphylo].
#' @export
plot.aphylo_estimates <- function(
  x,
  y = NULL,
  which.tree = 1L,
  ids        = list(1:Ntip(x)[which.tree]),
  loo        = TRUE,
  nsamples   = 1L,
  ncores     = 1L,
  centiles   = c(.025, .5, .975),
  cl         = NULL,
  ...
) {
  
  if (inherits(x$dat, "multiAphylo")) {
    if (!(which.tree %in% seq_len(Ntrees(x))) | length(which.tree) > 1L)
      stop("`which.tree` out of range.", call. = FALSE)
  }
  
  # Computing predictions
  pred <- stats::predict(
    object     = x,
    which.tree = which.tree,
    ids        = ids,
    loo        = loo,
    nsamples   = nsamples,
    centiles   = centiles,
    ncores     = ncores,
    cl         = cl
  )
  
  pred <- if (is.multiAphylo(x$dat))
    pred[[1]][1:Ntip(x)[which.tree], ,drop = FALSE]
  else
    pred[1:Ntip(x)[which.tree], ,drop = FALSE]
  
  # Adding the Pred. predix to the columns.
  colnames(pred) <- paste("Pred.", colnames(pred))
  
  if (is.multiAphylo(x$dat))
    x$dat <- x$dat[[which.tree]]
  
  x$dat$tip.annotation <- cbind(
    x$dat$tip.annotation,
    predicted = pred
  )
  
  if (nsamples > 1L)
    plot(x$dat, as_ci = 2L:4L, ...)
  else
    plot(x$dat, ...)
  
  invisible(x$dat)
}
