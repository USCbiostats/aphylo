#' Leave-one-out Cross Validation
#' 
#' This implements Leave-one-out cross-validation (LOO-CV) for trees of class
#' [aphylo] and [multiAphylo].
#' 
#' @param model As passed to [aphylo_mcmc].
#' @param ... Further arguments passed to the method.
#' @return An object of class `aphylo_cv` with the following components:
#' - `pred_out` Out of sample prediction.
#' - `expected` Expected annotations
#' - `call` The call
#' - `ids` Integer vector with the ids of the leafs used in the loo process.
#' 
#' @details For each observation in the dataset (either a single gene if of 
#' class [aphylo], or an entire tree if of class [multiAphylo]), we restimate
#' the model removing the observation and use the parameter estimates to make
#' a prediction on it. The prediction is done using the function [predict.aphylo_estimates]
#' with argument `loo = TRUE`.
#'  
#' @export
#' @examples 
#' # It takes about two minutes to run this example
#' \donttest{
#' 
#'   set.seed(123)
#'   atrees <- rmultiAphylo(10, 10, P = 1)
#' 
#'   cv_multi  <- aphylo_cv(atrees ~ mu_d + mu_s + Pi)
#'   cv_single <- aphylo_cv(atrees[[1]] ~ mu_d + mu_s + Pi)
#'   
#' }
aphylo_cv <- function(...) UseMethod("aphylo_cv")

#' @export
#' @rdname aphylo_cv
aphylo_cv.formula <- function(model, ...) {
  
  # First run of the model
  ans0    <- aphylo_mcmc(model, ...)
  ntrees  <- Ntrees(ans0)
  if (ntrees == 1) {
    has_ann <- which(rowSums(ans0$dat$tip.annotation == 9) < Nann(ans0$dat))
    nhas    <- length(has_ann)
  } else {
    nhas    <- ntrees
    has_ann <- 1L:nhas
  }
  
  # Model
  m <- as.formula(model)
  environment(m) <- environment()
  m[[2]] <- bquote(tree1)
  
  cat(sprintf("%s\nLeave-one-out cross validation of aphylo model with %i cases\n",
                 paste0(rep("-", 80L), collapse=""), nhas))
  pcents <- floor((1L:nhas)/nhas*100)
  
  # Output matrix
  if (ntrees == 1L) {
    
    pred <- matrix(
      9L,
      ncol = Nann(ans0$dat),
      nrow = Nnode(ans0$dat, internal.only = FALSE),
      # Proper row and column names
      dimnames = list(
        c(rownames(ans0$dat$tip.annotation), rownames(ans0$dat$node.annotation)),
        colnames(ans0$dat$node.annotation)
      )
    )
    
  } else {
    
    pred <- vector("list", ntrees)
    pred <- lapply(ans0$dat, function(tree.) {
      
      matrix(
        9L,
        ncol = Nann(tree.),
        nrow = Nnode(tree., internal.only = FALSE),
        # Proper row and column names
        dimnames = list(
          c(rownames(tree.$tip.annotation), rownames(tree.$node.annotation)),
          colnames(tree.$node.annotation)
        )
      )
      
    })
    
  }
  
  # Figuring out the iteration sequence
  iterseq <- seq_len(nhas)
  
  for (i in iterseq) {
    
    # Getting alternative model
    tree1 <- ans0$dat
    if (ntrees == 1L) {
      tree1[has_ann[i],] <- NA
    } else {
      tree1 <- tree1[-i]
    }
    
    ans1 <- suppressWarnings(suppressMessages(aphylo_mcmc(m, ...)))
    
    if (ntrees == 1) 
      pred[has_ann[i],] <- predict.aphylo_estimates(ans1)[has_ann[i],,drop=FALSE]
    else 
      pred[[i]] <- predict.aphylo_estimates(ans1, newdata = ans0$dat[[i]])
    
    
    # Communicating status
    if (interactive())
      cat(sprintf("\r %i of %i (% 3i%%) done...%s", i, nhas, pcents[i], c("\\", "/")[1 + i %% 2]))
    else
      message(sprintf("% 3i done...", has_ann[i]), appendLF = FALSE)
    
  }
  
  expected <- if (ntrees == 1) 
    with(ans0$dat, rbind(tip.annotation, node.annotation))
  else 
    lapply(ans0$dat, function(d) rbind(d$tip.annotation, d$node.annotation))
  
  
  structure(
    list(
      pred_out  = pred,
      expected  = expected,
      call      = sys.call(),
      ids       = iterseq,
      estimates = ans0,
      auc       = if (ntrees == 1) {
        auc(pred, expected)
        } else Map(auc, pred = pred, labels = expected),
      pscore    = if (ntrees == 1) {
        prediction_score(pred, expected)
        } else Map(prediction_score, x = pred, expected = expected),
      ntrees    = Ntrees(ans0)
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

