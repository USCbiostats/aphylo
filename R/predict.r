
#' @rdname aphylo_estimates-class
#' @return In the case of the `predict` method, a `P` column numeric matrix
#' with values between \eqn{[0,1]} (probabilities).
#' @export
predict.aphylo_estimates <- function(object, ...) {
  
  # Running prediction function
  pred <- with(
    object,
    predict_pre_order(
      atree = dat,
      psi   = par[c("psi0", "psi1")],
      mu    = par[c("mu0", "mu1")],
      eta   = par[c("eta0", "eta1")],
      Pi    = par["Pi"]
    )
    )
  
  # Adding names
  dimnames(pred) <- list(
    with(object$dat$tree, c(tip.label, node.label)),
    colnames(object$dat$tip.annotation))
  
  pred
}

#' @rdname aphylo_estimates-class
#' @param expected Integer vector of length \eqn{n}. Expected values (either 0 or 1).
#' @param alpha Numeric scalar. Prior belief of the parameter of the bernoulli distribution
#' used to compute the random imputation score.
#' @param W A square matrix. Must have as many rows as genes in `expected`.
#' @export
#' @details In the case of `prediction_score`, `...` are passed to
#' `predict.aphylo_estimates`.
#' 
#' @examples 
#' # Example with prediction_score ---------------------------------------------
#' set.seed(1312)
#' ap  <- sim_annotated_tree(10, P = 1, Pi=.2, mu=c(.05,.02))
#' ans <- aphylo_mcmc(rep(.05, 7), ap, control = list(nbatch=1e4, thin=100),
#'                    priors = function(x) dbeta(x, 1, 30))
#'                    
#' pr <- prediction_score(ans)
#' with(pr, cbind(Expected = expected, Predicted = predicted))
prediction_score <- function(
  x,
  expected = NULL,
  alpha    = 0.5,
  W        = NULL,
  ...) {
  
  # Finding relevant ids
  if (!length(expected)) 
    expected <- with(x$dat, rbind(tip.annotation, node.annotation))
  else {
    test <- all(dim(expected) == dim(with(x$dat, rbind(tip.annotation, node.annotation))))
    if (!test) 
      stop(
        "`expected` must be a matrix of size ",
        length(with(x$dat$tree, c(tip.label, node.label))), "x",
        ncol(x$dat$node.annotation), call. = FALSE)
  }

  # We will only focuse on those that we can actually asses
  ids <- which(apply(expected, 1L, function(x) all(x != 9L)))
  
  # And furthermore, only on the leafs
  ids <- intersect(ids, 1L:nrow(x$dat$tip.annotation))

  # Prediction
  pred <- predict.aphylo_estimates(x, ...)[ids,,drop=FALSE]
  
  # Inverse of Geodesic distances
  if (!length(W)) {
    G     <- approx_geodesic(x$dat$tree$edge - 1L, undirected = TRUE)[ids,ids]
    G_inv <- 1/(G + diag(nrow(G)))
  } else {
    G_inv <- W
    if (!all(dim(W) == rep(length(ids), 2)))
      stop(sprintf("-W- must have be of dimmension dim(W) == c(%i, %1$i)", length(ids)))
  }
  
  
  # Observed score
  if (!length(expected))
    expected <- x$dat$annotations[ids, ]
  
  obs <- sqrt(rowSums((pred - expected[ids,,drop = FALSE])^2))
  obs <- t(obs) %*% G_inv %*% obs
  
  # Best case
  best <- 0
  
  # Worst case
  worse <- sum(G_inv)*ncol(pred)
  
  # Random case
  rand  <- prediction_score_rand(expected[ids,,drop=FALSE], G_inv, alpha)
  
  structure(
    list(
      obs       = obs,
      worse     = worse,
      predicted = pred,
      expected  = expected[ids, ,drop=FALSE],
      random    = rand,
      alpha     = alpha
    ), class = "aphylo_prediction_score"
  )
  
}

predict_random <- function(P, A, G_inv) {
  n <- nrow(G_inv)
  sapply(1:10000, function(x) {
    A_hat <- matrix(sample(c(0,1), P*n, TRUE), ncol = P)
    obs   <- sqrt(rowSums((A - A_hat)^2))
    t(obs) %*% G_inv %*% obs
  })
}

#' @export
#' @rdname aphylo_estimates-class
print.aphylo_prediction_score <- function(x, ...) {
  cat("PREDICTION SCORE: ANNOTATED PHYLOGENETIC TREE\n")
  with(x, cat(
    sprintf("Observed : %-.2f ", obs/worse),
    sprintf("Random   : %-.2f ", random/worse),
    paste0(rep("-", getOption("width")), collapse=""),
    "Values scaled to range between 0 and 1, 0 being best.",
    sep ="\n"
  ))
  invisible(x)
}


#' @export
#' @param y Ignored.
#' @param main Passed to `title`.
#' @param which.fun Integer vector. Which function to plot.
#' @param include.labels Logical scalar. When `TRUE`, draws nice labels
#' at each slice which by default are specified as the rownames of `x$expected`.
#' This is mostly useful when the number of predictions is small.
#' @param labels.col Character scalar. Color of the labels.
#' @param main.colorkey Character scalar. Title of the colorkey (optional).
#' @rdname aphylo_estimates-class  
#' @details If `include.labels = NULL` and `ncol(x$expected) > 40`,
#' then `include.labels=FALSE` by default.
#' 
plot.aphylo_prediction_score <- function(
  x,
  y=NULL, 
  main = "Predicted v/s\nExpected Values",
  main.colorkey = "Probability of Functional Annotation",
  which.fun = seq_len(ncol(x$expected)),
  include.labels = NULL,
  labels.col = "black",
  ...) {
  
  k <- length(which.fun)
  y <- rep(1L, nrow(x$expected))
  
  # Should we draw the labels?
  if (!length(include.labels)) {
    if (nrow(x$expected) > 40) include.labels <- FALSE
    else include.labels <- TRUE
  }
    
  
  oldpar <- graphics::par(mfrow=c(1, k), mar=c(3,0,3,0))
  on.exit(graphics::par(oldpar))
  
  for (i in 1:k) {
    
    # Sorting accordingly to predicted
    ord <- 1L:length(x$predicted[,i]) # order(x$predicted[,i])

    
    # Outer polygon
    piechart(
      y, border="transparent", col = "transparent", lwd=2,
      radius = 1.5,
      doughnut = .5, skip.plot.slices = TRUE
      )
    
    graphics::polygon(polygons::circle(0,0,1.5), border="gray", lwd = 1.5)
    graphics::polygon(polygons::circle(0,0,0.5), border="gray", lwd = 1.5)

    # Function to color the absence/presence of function
    blue <- function(x) {
      ans <- polygons::colorRamp2(RColorBrewer::brewer.pal(7, "RdBu"))(x)
      grDevices::rgb(ans, alpha = 200, maxColorValue = 255)
    }
    
    # Outer pie
    polygons::piechart(
      y,
      radius    = 1,
      doughnut  = .755,
      add       = TRUE,
      col       = blue(x$predicted[ord,i]),
      border    = blue(x$predicted[ord,i]), 
      lwd       = 1.5,
      slice.off = abs(x$predicted[ord, i] - x$expected[ord, i])/2
    )
    
    # Inner pie
    polygons::piechart(
      y,
      doughnut  = 0.5,
      radius    = .745,
      add       = TRUE,
      col       = blue(x$expected[ord,i]),
      border    = blue(x$expected[ord,i]),
      lwd       = 1.5,
      slice.off = abs(x$predicted[ord, i] - x$expected[ord, i])/2
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
        labels    = rownames(x$expected)[ord],
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
    
    # graphics::polygon(circle(0, 0, r=.60), border = "black", lwd=1)
    graphics::text(0, 0, label=colnames(x$expected)[i], font=2)
  }
  
  # Drawing color key
  oldmar <- graphics::par(mar = rep(0, 4))
  graphics::par(mfrow=c(1,1))
  
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
  
