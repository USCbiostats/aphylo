rm(list = ls())

load("playground/simulations/data_and_functions.rda")
library(aphylo)

set.seed(1223)
ans_MAP <- vector("list", nsim)
for (i in 1:nsim) {
  
  # MAP estimators
  ans_MAP[[i]] <- mle_lite(dat_obs[[i]], FALSE, priors = function(params) dbeta(params, 1, 20))
  
  # Printing on screen (and saving)
  if (!(i %% 10)) {
    message(sprintf("Simulation %04i/%04i (%6.2f %%) complete.", i, nsim, i/nsim*100))
    
    # Just in case we need to restart this... so we can continue where we were
    curseed_MAP <- .Random.seed
    
    # Storing all objects
    save(curseed_MAP, ans_MAP, i, file = "playground/simulations/map_estimates.rda")
  } else message(".", appendLF = FALSE)
  
}

# Terminating
message(sprintf("Simulation %04i/%04i (%6.2f %%) complete.", i, nsim, i/nsim*100))

curseed_MAP <- .Random.seed
save(curseed_MAP, ans_MAP, i, file = "playground/simulations/map_estimates.rda")
