
#' @rdname mle
#' @param what Either a character scalar or an integer vector. If a character,
#' then it can be either \code{"missings"}, \code{"leafs"}, or \code{"all"}. If an integer vector,
#' then these must be values between \eqn{[0, n - 1]} (node ids).
#' @return In the case of the \code{predict} method, a two-column numeric matrix
#' with values between \eqn{[0,1]} (probabilities).
#' @export
predict.phylo_mle <- function(object, what = c("missings", "all"), ...) {
  
  # Parameters
  n <- nrow(object$dat$annotations)
  
  # Checking the default
  if (length(what) == 2 && all(what == c("missings", "all")))
    what <- "missings"
  
  # Checking what to predict
  if (length(what) >= 1 && inherits(what, "integer")) {
    ran <- range(what)
    
    # Out of range
    test <- which(ran < 0 | ran >= n)
    if (length(test))
      stop("Ids in -what- out of range:\n", paste(what[test], collapse=", "), ".")
    
    ids <- what
  } else if (length(what) == 1 && what == "missings") {
    ids <- with(object$dat, which(
      noffspring == 0 & apply(annotations, 1, function(a) any(a == 9)*1L) > 0
    ))
    
    if (!length(ids))
      stop("No missing nodes to predict.")
    
    # Adjusting indices
    ids <- ids - 1L
    
  } else if (length(what) == 1 && what == "all") {
    ids <- 0L:(n-1L)
  } else if (length(what) == 1 && what == "leafs") {
    ids <- which(object$dat$noffspring == 0) - 1L
  }else 
    stop("Undefined method for -what- equal to: ", what)
  
  # Running prediction function
  pred <- with(object, 
               predict_funs(
                 ids         = ids,
                 edges       = dat$edges,
                 annotations = dat$annotations,
                 offspring   = dat$offspring,
                 noffspring  = dat$noffspring,
                 psi         = par[1:2],
                 mu          = par[3:4],
                 Pi          = c(1 - par[5], par[5])
               )
  )
  
  # Adding names
  dimnames(pred) <- list(ids, colnames(object$dat$annotations))
  
  pred
}

#' @rdname mle
#' @param expected Integer vector of length \eqn{n}. Expected values (either 0 or 1).
#' @export
#' @details In the case of \code{prediction_score}, \code{...} are passed to
#' \code{predict.phylo_mle}.
prediction_score <- function(x, expected = NULL, ...) {
  
  # Finding relevant ids
  if (!length(expected)) 
    expected <- x$dat$annotations
  else {
    test <- all(dim(expected) == dim(x$dat$annotations))
    if (!test) 
      stop("-expected- must have the same dimmension as -x$dat$annotations-.")
  }

  # We will only focuse on those that we can actually asses
  ids <- which(apply(expected, 1L, function(x) all(x != 9L)))
  
  # And furthermore, only on the leafs
  ids <- intersect(ids, which(x$dat$noffspring == 0))

  # Prediction
  pred <- predict.phylo_mle(x, what = ids - 1L, ...)
  
  # Inverse of Geodesic distances
  G     <- approx_geodesic(x$dat$edges, undirected = TRUE)[ids,ids]
  G_inv <- 1/(G + 1)
  
  # Observed score
  if (!length(expected))
    expected <- x$dat$annotations[ids, ]
  
  obs <- sqrt(rowSums((pred - expected[ids,,drop = FALSE])^2))
  obs <- t(obs) %*% G_inv %*% obs
  
  # Best case
  best <- 0
  
  # Worst case
  worse <- sum(G_inv)*ncol(pred)
  
  structure(
    list(
      obs       = obs,
      worse     = worse,
      predicted = pred,
      expected  = expected[ids, ,drop=FALSE]
    ), class = "aphylo_prediction_score"
  )
  
}

#' @export
#' @rdname mle
print.aphylo_prediction_score <- function(x, ...) {
  cat("PREDICTION SCORE: ANNOTATED PHYLOGENETIC TREE\n")
  with(x, cat(sprintf("Observed: %.2f\nBest: 0.00\nWorse: %.2f\nRelative (obs/worse):%.2f\n",
                      obs, worse, 1 - obs/worse)
              )
       )
  invisible(x)
}
  
  