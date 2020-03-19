#' Calculate prediction score (quality of prediction)
#' 
#' @param x An object of class [aphylo_estimates] or a numeric matrix.
#' @param expected Integer vector of length \eqn{n}. Expected values (either 0 or 1).
#' @param alpha Numeric scalar. Prior belief of the parameter of the bernoulli distribution
#' used to compute the random imputation score.
#' @param W A square matrix. Must have as many rows as genes in `expected`.
#' @param ... Further arguments passed to [predict.aphylo_estimates]
#' @export
#' @details In the case of `prediction_score`, `...` are passed to
#' `predict.aphylo_estimates`.
#' 
#' @examples 
#' # Example with prediction_score ---------------------------------------------
#' set.seed(11552)
#' ap  <- raphylo(50, P = 1, Pi=0, mu_d=c(.8,.2), mu_s = c(0,0.25), psi = c(0,0))
#' ans <- aphylo_mcmc(
#'   ap ~ mu_d + mu_s + psi + Pi,
#'   control = list(nsteps=5e3, thin=20, burnin = 1000),
#'   priors = bprior(c(9,9,1,1,1,1,1), c(1,1,9,9,9,9,9))
#'   )
#'                    
#' pr <- prediction_score(ans)
#' with(pr, cbind(Expected = expected, Predicted = predicted))
prediction_score <- function(x, expected, alpha = NULL, W = NULL, ...)
  UseMethod("prediction_score")

#' @export
#' @rdname prediction_score
prediction_score.default <- function(x, expected, alpha = NULL, W = NULL, ...) {
  
  # Checking dimensions
  if (length(x) != length(expected))
    stop("`x` and `expected` differ in length. These must match.", call.=FALSE)
  
  # Counting complete cases
  N <- 1:nrow(expected)
  ids <- expected[cbind(N, max.col(expected))]
  ids <- which(ids %in% c(1L, 0L) & x[cbind(N, max.col(x))] != 9.0)
  
  if (length(ids) != nrow(expected)) {
    expected <- expected[ids, , drop=FALSE]
    x        <- x[ids, , drop=FALSE]
  }
  
  # Computing the expected value of 1s. We will use this to compute the random
  # score.
  if (!length(alpha))
    alpha <- mean(expected)
  
  if (is.null(W))
    W <- diag(length(ids))
  else
    W <- W[ids, ids, drop=FALSE]
  
  obs <- rowSums(x - expected)
  obs <- t(obs) %*% W %*% obs
  
  # Best case
  best <- 0
  
  # Worst case
  worse <- sum(W)*ncol(x)
  
  # Random case
  rand  <- prediction_score_rand(expected, W, alpha)
  
  # Computing AUCs
  
  
  structure(
    list(
      obs       = obs,
      worse     = worse,
      predicted = x,
      expected  = expected,
      random    = rand,
      alpha     = alpha,
      auc       = auc(x, expected),
      obs.ids   = NULL,
      leaf.ids  = NULL,
      tree      = NULL
    ), class = "aphylo_prediction_score"
  )
  
}

#' @export
#' @template loo
#' @rdname prediction_score
#' @details In the case of the method for aphylo estimates, the function takes as
#' a reference using alpha equal to the proportion of observed tip annotations that
#' are equal to 1, this is:
#' 
#' ```
#' mean(x$dat$tip.annotation[x$dat$tip.annotation != 9L], na.rm = TRUE)
#' ```
prediction_score.aphylo_estimates <- function(
  x,
  expected = NULL,
  alpha    = NULL,
  W        = NULL,
  loo      = TRUE,
  ...
  ) {
  
  # If null, then we take the true (observed) proportion of ones to compare
  if (any(class(x$dat) %in% c("multiAphylo"))) {
    
    # Setting up for prediction
    x_tmp <- x
    res   <- vector("list", Ntrees(x_tmp))
    
    if (is.null(expected)) expected <- res
    if (is.null(W))        W        <- res
    if (length(alpha) == 1L)  alpha <- rep(alpha, Ntrees(x))
      
    for (i in seq_along(res)) {
      
      # Preparing data
      x_tmp$dat <- x$dat[[i]]
      res[[i]] <- prediction_score.aphylo_estimates(
        x_tmp,
        expected = expected[[i]],
        alpha    = alpha[i],
        W        = W[[i]],
        loo      = loo,
        ...
        )
    }
    
    return(res)
    
  }
  
  if (is.null(alpha))
    alpha <- mean(x$dat$tip.annotation[x$dat$tip.annotation != 9L], na.rm = TRUE)
  
  # Finding relevant ids
  if (is.null(expected)) {
    expected <- with(x$dat, rbind(tip.annotation, node.annotation))
    dimnames(expected) <- with(x$dat, list(c(tree$tip.label, tree$node.label), colnames(tip.annotation)))
  } else {
    test <- all(dim(expected) == dim(with(x$dat, rbind(tip.annotation, node.annotation))))
    if (!test) 
      stop(
        "`expected` must be a matrix of size ",
        length(with(x$dat$tree, c(tip.label, node.label))), "x",
        Nann(x$dat), call. = FALSE)
  }

  # We will only focuse on those that we can actually asses
  ids <- which(apply(expected, 1L, function(x) all(x != 9L)))
  
  # And furthermore, only on the leafs
  ids <- intersect(ids, 1L:Ntip(x$dat))

  # Prediction
  pred <- predict.aphylo_estimates(x, ids = list(ids), loo = loo,...)
  
  # Inverse of Geodesic distances
  if (!length(W)) {
    # G     <- approx_geodesic(x$dat$tree$edge - 1L, undirected = TRUE)[ids,ids]
    G_inv <- diag(length(ids))
  } else {
    G_inv <- W
    if (!all(dim(W) == rep(length(ids), 2)))
      stop(sprintf("-W- must have be of dimmension dim(W) == c(%i, %1$i)", length(ids)))
  }
  
  ans <- prediction_score(
    x        = pred[ids,,drop=FALSE],
    expected = expected[ids,,drop = FALSE],
    alpha    = alpha,
    W        = G_inv
  )
  
  # Adding missing info
  ans$obs.ids   <- c(x$dat$tree$tip.label,x$dat$tree$node.label)[ids]
  ans$tree      <- x$dat$tree
  
  ans$predicted <- pred
  ans$expected  <- expected
  
  ans
  
}

#' Calculates the ramdon prediction score by simulationss
#' @param P Number of functions
#' @param A Observed annotations
#' @param G_inv Weighting matrix
#' @noRd
#' 
predict_random <- function(P, A, G_inv, alpha, R = 1e4L) {
  n <- nrow(G_inv)
  sapply(1:R, function(x) {
    
    A_hat <- matrix(
      data = sample(
        x       = c(0, 1),
        size    = P * n,
        replace = TRUE,
        prob    = c(1 - alpha, alpha)
        ),
      ncol = P
      )
    
    obs   <- sqrt(rowSums((A - A_hat)^2))
    t(obs) %*% G_inv %*% obs
  })
}

#' @export
#' @rdname prediction_score
print.aphylo_prediction_score <- function(x, ...) {
  cat("PREDICTION SCORE: ANNOTATED PHYLOGENETIC TREE\n")
  with(x, cat(
    sprintf("Observed : %-.2f ", obs/worse),
    sprintf("Random   : %-.2f ", random/worse),
    sprintf("AUC      : %-.2f ", auc$auc),
    paste0(rep("-", getOption("width")), collapse=""),
    "Values scaled to range between 0 and 1, 0 being best.",
    sep ="\n"
  ))
  invisible(x)
}


# Function to color the absence/presence of function
blue <- function(x) {
  x[is.na(x)] <- 9
  ans <- polygons::colorRamp2(.aphyloColors)(x)
  ans <- grDevices::rgb(ans, alpha = 200, maxColorValue = 255)
  ifelse(x == 9, "black", ans)
}

#' Visualize predictions
#' 
#' @export
#' @param x An object of class `aphylo_prediction_score`.
#' @param ... Ignored
#' @param y Ignored.
#' @param main Passed to `title`.
#' @param which.fun Integer vector. Which function to plot.
#' @param include.labels Logical scalar. When `TRUE`, draws nice labels
#' at each slice which by default are specified as the rownames of `x$expected`.
#' This is mostly useful when the number of predictions is small.
#' @param labels.col Character scalar. Color of the labels.
#' @param main.colorkey Character scalar. Title of the colorkey (optional).
#' @param leafs_only Logical. When `TRUE` (default) only plots the leaf nodes.
#' @details
#' 
#' If `include.labels = NULL` and `ncol(x$expected) > 40`,
#' then `include.labels=FALSE` by default.
#' @aliases plot-prediction
plot.aphylo_prediction_score <- function(
  x,
  y              = NULL, 
  main           = "Prediction Accuracy: Observed versus predicted values",
  main.colorkey  = "Probability of Functional Annotation",
  which.fun      = seq_len(ncol(x$expected)),
  include.labels = NULL,
  labels.col     = "black",
  leafs_only     = TRUE,
  ...
  ) {
  
  if (is.null(x$tree))
    stop("This method is only available for trees.", call. = FALSE)
  
  # Should we plot only the leafs?
  if (leafs_only) {
    idx <- x$tree$tip.label
  } else
    idx <- 1:nrow(x$predicted)
  
  predicted <- x$predicted[idx, , drop=FALSE]
  expected  <- x$expected[idx, , drop=FALSE]
  
  k <- length(which.fun)
  y <- rep(1L, nrow(predicted))
  
  # Should we draw the labels?
  if (!length(include.labels)) {
    
    if (nrow(predicted) > 40) include.labels <- FALSE
    else include.labels <- TRUE
    
  }
  
  # Getting the order
  grDevices::pdf(file = NULL);ape::plot.phylo(x$tree, plot=FALSE);grDevices::dev.off()
  plot_pars <- utils::getFromNamespace(".PlotPhyloEnv", "ape")
  ord  <- order(predicted[,1])#order(plot_pars$last_plot.phylo$yy[idx])
  
  oldpar <- graphics::par(mar=c(3,0,3,0))
  on.exit(graphics::par(oldpar))
  
  for (i in 1:k) {
    
    # Sorting accordingly to predicted
    # ord <- 1L:length(predicted[,i]) 

    # Outer polygon
    piechart(
      y, border="transparent", col = "transparent", lwd=2,
      radius = 1.5,
      doughnut = .5, skip.plot.slices = TRUE
      )
    
    graphics::polygon(polygons::circle(0,0,1.5), border="gray", lwd = 1.5, col = "lightgray")
    graphics::polygon(polygons::circle(0,0,0.5), border="gray", lwd = 1.5, col="white")
    
    # Outer pie
    opie <- polygons::piechart(
      y,
      radius    = 1,
      doughnut  = .755,
      add       = TRUE,
      col       = blue(predicted[ord,i]),
      border    = blue(predicted[ord,i]), 
      lwd       = .5,
      slice.off = ifelse(
        expected[ord, i] == 9L,.25,
        abs(predicted[ord, i] - expected[ord, i])/2
      )
    )
    
    # Inner pie
    ipie <- polygons::piechart(
      y,
      doughnut  = 0.5,
      radius    = .745,
      add       = TRUE,
      col       = blue(expected[ord,i]),
      border    = blue(expected[ord,i]),
      lwd       = .5,
      density   = ifelse(expected[ord,i] == 9L, 10, NA),
      slice.off = ifelse(
        expected[ord, i] == 9L,.25,
        abs(predicted[ord, i] - expected[ord, i])/2
      )
    )
    
    # Labels
    if (include.labels) {
      deg <- 1:length(y)
      deg <- c(deg[1], deg[-1] + deg[-length(y)])/length(y)/2*360
      polygons::piechart(
        y,
        radius    = .5,
        add       = TRUE,
        border    = NA,
        labels    = rownames(expected)[ord],
        text.args = list(
          srt  = ifelse(deg > 270, deg,
                        ifelse(deg > 90, deg + 180, deg)),
          col  = labels.col, #c("white", "darkgray"),
          cex  = .7 - (1 - 1/k)*.5,
          xpd  = TRUE
        ),
        tick.len  = 0,
        segments.args = list(col="transparent"),
        skip.plot.slices = TRUE
      )
    }
    
    # Extra annotations
    graphics::segments(
      x0 = 0.1, y0 = 1.6, x1 = 1.75, y1 = 1.6, col="black", 
      lty=2, lwd=2
    )
    graphics::segments(x0 = 0, y0 = 1.5, x1 = .1, y1 = 1.6, lwd=2)
    graphics::text(1.75,1.6, labels = "Perfect miss", pos = 3)
    
    graphics::segments(
      x0 = 0.1, y0 = 0.6, x1 = 1.75, y1 = 0.6, col="black", 
      lty=2, lwd=2)
    graphics::segments(x0 = 0, y0 = .5, x1 = .1, y1 = .6, lwd=2)
    graphics::text(1.75,0.6, labels = "Perfect\nprediction", pos = 3)
    
    # Adding more notes
    slice2annotate <- which.min(opie$textcoords[,1])
    
    opie <- colMeans(opie$slices[[slice2annotate]])
    ipie <- colMeans(ipie$slices[[slice2annotate]])
    
    graphics::text(-1.76, .7, labels = "Observed\nannotation", pos=3)
    graphics::segments(-1.76, .7, ipie[1], ipie[2], lty=2, lwd=2)
    graphics::text(-1.76, -.7, labels = "Predicted\nannotation", pos=1)
    graphics::segments(-1.76, -.7, opie[1], opie[2], lty=2, lwd=2)
    
    graphics::text(0, 0, label=colnames(expected)[i], font=2)
    
    # Drawing color key
    oldmar <- graphics::par(mar = rep(0, 4), new = FALSE, xpd=NA)
    # graphics::par(mfrow=c(1,1), xpd=NA)
    
    polygons::colorkey(
      x0 = .10, y0=0, x1=.90, y1=.1, 
      label.from = 'No function',
      label.to = "Function",
      cols = blue(seq(0,1,length.out = 3)), 
      tick.range = c(0,1),
      tick.marks = c(0,.25,.5,.75,1),
      nlevels = 200,
      main = main.colorkey
    )
    graphics::par(oldmar)
    
    graphics::title(
      main= main, font.main=1
    )
  }

}
  