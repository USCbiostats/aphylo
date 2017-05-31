rm(list = ls())
library(aphylo)
load("playground/simulations/data_and_functions.rda")

set.seed(1223)
ans_MLE <- vector("list", nsim)
for (i in 1:nsim) {
  
  # MLE estimator
  ans_MLE[[i]] <- mle_lite(dat_obs[[i]], FALSE)
  
  # Printing on screen (and saving)
  if (!(i %% 10)) {
    message(sprintf("Simulation %04i/%04i (%6.2f %%) complete.", i, nsim, i/nsim*100))
    
    # Just in case we need to restart this... so we can continue where we were
    curseed_MLE <- .Random.seed
    
    # Storing all objects
    save(curseed_MLE, ans_MLE, i, file = "playground/simulations/mle_estimates.rda")
  } else message(".", appendLF = FALSE)
  
}

# Terminating
message(sprintf("Simulation %04i/%04i (%6.2f %%) complete.", i, nsim, i/nsim*100))

curseed_MLE <- .Random.seed
save(curseed_MLE, ans_MLE, i, file = "playground/simulations/mle_estimates.rda")
