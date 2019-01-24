loo <- function(model, ...) {
  
  ans0 <- aphylo_mcmc(model, ...)
  nt   <- Ntip(ans0$dat)
  
  # Model
  m <- as.formula(model, env = sys.frame())
  m[[2]] <- bquote(tree1)
  
  pred <- vector("list", nt)
  for (i in 1L:nt) {
    
    # Getting alternative model
    tree1 <- ans0$dat
    tree1[i,] <- NA
    
    ans1 <- suppressWarnings(suppressMessages(aphylo_mcmc(m, ...)))
    pred[[i]] <- predict(ans1)[i,,drop=FALSE]
    
    message(sprintf("% 3i done...", i))
    
  }
  
  list(
    loo  = do.call(rbind, pred),
    mcmc = ans0,
    call = sys.call()
  )
    
  
}

x <- sim_annotated_tree(5)
ans <- loo(x ~ psi + mu + Pi, priors = bprior())
