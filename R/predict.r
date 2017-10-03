
#' @rdname aphylo_estimates-class
#' @param what Either a character scalar or an integer vector. If a character,
#' then it can be either \code{"missings"}, \code{"leafs"}, or \code{"all"}. If an integer vector,
#' then these must be values between \eqn{[0, n - 1]} (node ids).
#' @return In the case of the \code{predict} method, a two-column numeric matrix
#' with values between \eqn{[0,1]} (probabilities).
#' @export
predict.aphylo_estimates <- function(object, what = c("missings", "all"), ...) {
  
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
      isleaf(edges) & apply(annotations, 1, function(a) any(a == 9)*1L) > 0
    ))
    
    if (!length(ids))
      stop("No missing nodes to predict.")
    
    # Adjusting indices
    ids <- ids - 1L
    
  } else if (length(what) == 1 && what == "all") {
    ids <- 0L:(n-1L)
  } else if (length(what) == 1 && what == "leafs") {
    ids <- which(isleaf(object$dat$edges)) - 1L
  } else if (is.vector(what) & inherits(what, "character")) {
    
    # Fetching nodes using labels
    ids  <- match(what, rownames(object$dat$annotations))
    test <- which(is.na(ids))
    if (length(test))
      stop("Some elements of -what- are not present in the data: ",
           paste0(what[test], collapse=", "))
      
  } else 
    stop("Undefined method for -what- equal to: ", what)
  
  # Running prediction function
  pred <- with(object, 
               predict_funs(
                 ids         = ids,
                 edges       = dat$edges,
                 annotations = dat$annotations,
                 offspring   = attr(dat$edges, "offspring"),
                 psi         = par[1:2],
                 mu          = par[3:4],
                 Pi          = par[5]
               )
  )
  
  # Adding names
  dimnames(pred) <- list(
    rownames(object$dat$annotations)[ids+1L],
    colnames(object$dat$annotations))
  
  pred
}

#' @rdname aphylo_estimates-class
#' @param expected Integer vector of length \eqn{n}. Expected values (either 0 or 1).
#' @param alpha Numeric scalar. Prior belief of the parameter of the bernoulli distribution
#' used to compute the random imputation score.
#' @param W A square matrix. Must have as many rows as genes in \code{expected}.
#' @export
#' @details In the case of \code{prediction_score}, \code{...} are passed to
#' \code{predict.aphylo_estimates}.
#' 
#' @examples 
#' # Example with prediction_score ---------------------------------------------
#' set.seed(1312)
#' ap  <- sim_annotated_tree(10, P = 1, Pi=.2, mu=c(.05,.02))
#' ans <- aphylo_mcmc(rep(.05, 5), ap, control = list(nbatch=1e4, thin=100),
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
    expected <- x$dat$annotations
  else {
    test <- all(dim(expected) == dim(x$dat$annotations))
    if (!test) 
      stop("-expected- must have the same dimmension as -x$dat$annotations-.")
  }

  # We will only focuse on those that we can actually asses
  ids <- which(apply(expected, 1L, function(x) all(x != 9L)))
  
  # And furthermore, only on the leafs
  ids <- intersect(ids, which(isleaf(x$dat$edges)))

  # Prediction
  pred <- predict.aphylo_estimates(x, what = ids - 1L, ...)
  
  # Inverse of Geodesic distances
  if (!length(W)) {
    G     <- approx_geodesic(x$dat$edges, undirected = TRUE)[ids,ids]
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
    sprintf("Observed : %-.2f (%.2f)", obs/worse, obs),
    sprintf("Random   : %-.2f (%.2f)", random/worse, random),
    sprintf("Best     : 0.00 (0.00)"),
    sprintf("Worse    : 1.00 (%.2f)", worse),
    paste0(rep("-", getOption("width")), collapse=""),
    "Values between 0 and 1, 0 been best. Absolute scores in parenthesis.",
    sep ="\n"
  ))
  invisible(x)
}


#' @export
#' @param y Ignored.
#' @param main Passed to \code{title}.
#' @param which.fun Integer vector. Which function to plot.
#' @param include.labels Logical scalar. When \code{TRUE}, draws nice labels
#' at each slice which by default are specified as the rownames of \code{x$expected}.
#' This is mostly useful when the number of predictions is small.
#' @param labels.col Character scalar. Color of the labels.
#' @param main.colorkey Character scalar. Title of the colorkey (optional).
#' @rdname aphylo_estimates-class  
#' @details If \code{include.labels = NULL} and \code{ncol(x$expected) > 40},
#' then \code{include.labels=FALSE} by default.
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
    
  
  
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  graphics::par(mfrow=c(1, k), mar=c(6,0,3,0))
  
  for (i in 1:k) {
    
    # Sorting accordingly to predicted
    ord <- order(x$predicted[,i])

    
    # Outer polygon
    piechart(
      y, border="transparent", col = "transparent", lwd=2,
      radius = 1.5,
      doughnut = .5, skip.plot.slices = TRUE
      )
    
    graphics::polygon(circle(0,0,1.5), border="gray", lwd = 1.5)
    graphics::polygon(circle(0,0,0.5), border="gray", lwd = 1.5)

    # Function to color the absence/presence of function
    blue <- function(x) {
      ans <- grDevices::colorRamp(c("steelblue", "lightgray", "darkred"))(x)
      grDevices::rgb(ans[,1], ans[,2], ans[,3], 200, maxColorValue = 255)
    }
    
    # Outer pie
    piechart(
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
    piechart(
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
      piechart(
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
  graphics::par(mfrow=c(1,1), new=TRUE, mar=c(3,0,3,0), xpd=TRUE)
  graphics::plot.new()
  graphics::plot.window(c(0,1), c(0,1))
  
  colorkey(
    .10, 0, .90, .05, 
    label.from = 'No function',
    label.to = "Function",
    cols = grDevices::adjustcolor(c("steelblue", "lightgray", "darkred"), alpha.f = 200/255),
    tick.range = c(0,1),
    tick.marks = c(0,.25,.5,.75,1),
    nlevels = 200,
    main = main.colorkey
  )
  
  graphics::title(
    main= main, font.main=1
    )

}
  