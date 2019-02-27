#' Leave-one-out Cross Validation
#' @param model As passed to [aphylo_mcmc].
#' @param ... Further arguments passed to the method.
#' @return An object of class `aphylo_cv` with the following components:
#' - `pred_out` Out of sample prediction.
#' - `expected` Expected annotations
#' - `call` The call
#' - `ids` Integer vector with the ids of the leafs used in the loo process.
#' @export
aphylo_cv <- function(...) UseMethod("aphylo_cv")

#' @export
#' @rdname aphylo_cv
aphylo_cv.formula <- function(model, ...) {
  
  # First run of the model
  ans0    <- aphylo_mcmc(model, ...)
  nt      <- Ntip(ans0$dat)
  has_ann <- which(rowSums(ans0$dat$tip.annotation == 9) < Nann(ans0$dat))
  nhas    <- length(has_ann)
  
  # Model
  m <- as.formula(model)
  environment(m) <- environment()
  m[[2]] <- bquote(tree1)
  
  cat(sprintf("%s\nLeave-one-out cross validation of aphylo model with %i cases\n",
                 paste0(rep("-", 80L), collapse=""), nhas))
  pcents <- floor((1:nhas)/nhas*100)
  
  # Output matrix
  pred <- matrix(
    9L, ncol = Nann(ans0$dat), nrow = Nnode(ans0$dat, internal.only = FALSE),
    # Proper row and column names
    dimnames = list(
      c(rownames(ans0$dat$tip.annotation), rownames(ans0$dat$node.annotation)),
      colnames(ans0$dat$node.annotation)
    )
  )
  
  for (i in seq_along(has_ann)) {
    
    # Getting alternative model
    tree1 <- ans0$dat
    tree1[has_ann[i],] <- NA
    
    ans1 <- suppressWarnings(suppressMessages(aphylo_mcmc(m, ...)))
    pred[has_ann[i],] <- predict.aphylo_estimates(ans1)[has_ann[i],,drop=FALSE]
    
    # Communicating status
    if (interactive())
      cat(sprintf("\r %i of %i (% 3i%%) done...%s", i, nhas, pcents[i], c("\\", "/")[1 + i %% 2]))
    else
      message(sprintf("% 3i done...", has_ann[i]), appendLF = FALSE)
    
  }
  
  expected <- with(ans0$dat, rbind(tip.annotation, node.annotation))
  
  structure(
    list(
      pred_out  = pred,
      expected  = expected,
      call      = sys.call(),
      ids       = has_ann,
      estimates = ans0,
      auc       = auc(pred, expected),
      pscore    = prediction_score(pred, expected)
    ),
    class="aphylo_cv"
  )
  
}

#' @export
#' @param x An object of class `aphylo_auc`.
#' @param ... Further arguments passed to the method.
#' @rdname auc
print.aphylo_auc <- function(x, ...) {
  
  with(x, {
    cat(sprintf("Number of observations     : %i\n", n_used))
    cat(sprintf("Area Under The Curve (AUC) : %03.2f\n", auc))
    cat("Rates can be accessed via the $ operator.\n")
  })
  
  invisible(x)
  
}

#' @export
#' @param y Ignored.
#' @rdname auc
plot.aphylo_auc <- function(x, y=NULL, ...) {
  
  dots <- list(...)
  if (!length(dots$xlab)) dots$xlab <- "False Positive Rate"
  if (!length(dots$ylab)) dots$ylab <- "True Positive Rate"
  if (!length(dots$main)) dots$main <- "Receiver Operating Characteristic"
  if (!length(dots$type)) dots$type <- "l"
  
  dots$x <- x$fpr
  dots$y <- x$tpr
  
  do.call(graphics::plot, dots)
  graphics::abline(a=0, b=1, col = "gray", lty="dashed", lwd=1.5)
  
}