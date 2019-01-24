library(aphylo)

aphylo_cv <- function(...) UseMethod("aphylo_cv")

aphylo_cv.formula <- function(model, ...) {
  
  # First run of the model
  ans0    <- aphylo_mcmc(model, ...)
  nt      <- Ntip(ans0$dat)
  has_ann <- which(rowSums(ans0$dat$tip.annotation == 9) < Nann(ans0$dat))
  nhas    <- length(has_ann)
  
  # Model
  m <- as.formula(model, env = sys.frame())
  m[[2]] <- bquote(tree1)
  
  msg <- sprintf("%s\nLeave-one-out cross validation of aphylo model with %i cases",
                 paste0(rep("-", 80L), collapse=""), nhas)
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
    pred[has_ann[i],] <- predict(ans1)[has_ann[i],,drop=FALSE]
    
    # Communicating status
    if (interactive())
      cat(sprintf("\f%s\n %i of %i (% 3i%%) done...%s", msg, i, nhas, pcents[i], c("\\", "/")[1 + i %% 2]))
    else
      message(sprintf("% 3i done...", has_ann[i]), appendLF = FALSE)
    
  }
  
  structure(
    list(
      pred_out  = pred,
      expected  = with(ans0$dat, rbind(tip.annotation, node.annotation)),
      call      = sys.call(),
      ids       = has_ann
      ),
    class="aphylo_cv"
    )
    
  
}

# print.aphylo_cv

set.seed(123)
x <- sim_annotated_tree(20)
# x <- rdrop_annotations(x, .8, informative = TRUE)
ans <- aphylo_cv(x ~ psi + mu + Pi, priors = bprior())
# 
# "0        20        40        60        80       100"
# "|---------|---------|---------|---------|---------|"