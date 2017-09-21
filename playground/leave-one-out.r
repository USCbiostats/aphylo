rm(list = ls())

library(aphylo)

set.seed(1313)

dat<- sim_annotated_tree(n = 200, P = 1)

plot(dat)

ans <- aphylo_mcmc(rep(0.05, 5), dat, priors = function(param) dbeta(param, 1, 10),
                   control = list(nbatch=1e5, nchains=5))

(score <- prediction_score(ans))
plot(score)


dat

prune <- function(x, which) {
  leafs <- leafs(x)
  
}