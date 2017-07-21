rm(list = ls())

library(parallel)

# devtools::load_all()

set.seed(1122)

n          <- 2e3
parameters <- lapply(1:n, function(x) rbeta(5, 1, 10))# runif(5))
datasets   <- lapply(parameters, function(x) {
  aphylo::sim_annotated_tree(100, psi=x[1:2], mu=x[3:4], Pi = x[5])
})

try_phylo_mle <- function(...)
  tryCatch(phylo_mle(...), error = function(e) NA)

cl <- makeForkCluster(16)

# Old Version
system("git checkout master")
e0 <- clusterEvalQ(cl, {
  devtools::load_all()
})

ans_old <- parLapply(cl, datasets, try_phylo_mle)

# New version
system("git checkout multiple-start")
e0 <- clusterEvalQ(cl, {
  devtools::load_all()
  })

ans_new <- parLapply(cl, datasets, try_phylo_mle)

# Pi 0/1
system("git checkout mlePi01")
e0 <- clusterEvalQ(cl, {
  devtools::load_all()
  })

ans_Pi01 <- parLapply(cl, datasets, try_phylo_mle)

# Fixing pi
system("git checkout fixpi")
e0 <- clusterEvalQ(cl, {
  devtools::load_all()
})

ans_fixpi <- clusterMap(
  cl, try_phylo_mle,
  dat = datasets, 
  Pi  = lapply(1:n, function(x) parameters[[x]][5]))

stopCluster(cl)
# Tabulating results -----------------------------------------------------------

extractme <- function(x) {
  if (inherits(x, "phylo_mle")) x[["ll"]]
  else NA
}

mle_old   <- unlist(sapply(ans_old, extractme))
mle_new   <- unlist(sapply(ans_new, extractme))
mle_Pi01  <- unlist(sapply(ans_Pi01, extractme))
mle_fixpi <- unlist(sapply(ans_fixpi, extractme))

loglikes <- data.frame(
  old   = mle_old,
  new   = mle_new,
  Pi01  = mle_Pi01,
  fixpi = mle_fixpi,
  largest = apply(cbind(mle_old, mle_new, mle_Pi01, mle_fixpi), 1, function(x) {
    if (any(is.na(x))) return("-")
    paste(c("old", "No fixed", "Pi 0/1", "Pi = Pi*")[which(x == max(x))], collapse="-")
    
  }), check.names = FALSE
)

biascalc <- function(x, y) {
  if (inherits(x, "phylo_mle")) 
    x[["par"]] - y
  else rep(NA, 5)
}

biases <- rbind(
  data.frame(do.call(rbind, Map(biascalc, ans_old, parameters)), model = "MLE starting at {.05}^5"),
  data.frame(do.call(rbind, Map(biascalc, ans_new, parameters)), model = "MLE starting at {.05, .95}^5"),
  data.frame(do.call(rbind, Map(biascalc, ans_Pi01, parameters)), model = "MLE starting at {.05, .95}^4\n and comparing Pi 0/1"),
  data.frame(do.call(rbind, Map(biascalc, ans_fixpi, parameters)), model = "MLE starting at {.05, .95}^4\nandfixing Pi = Pi*")
)

# The new is larger than the old one,
table(loglikes$largest)

# But the size of the bias of the new one, is huge
oldpar <- par(no.readonly = TRUE)
par(mfrow = c(2, 2), oma = c(2,0,2,0))
invisible(by(biases, biases$model, function(x) {
  boxplot(x[,-6], ylim = c(-1, 1), main = unique(x$model))
  abline(h=0, lty="dashed", lwd=2)
}))
par(oldpar)
title(
  main="Bias distribution",
  sub = paste(
    paste("Each box represents", n, "simulations of annotated trees with 100 leafs."),
    "All parameters follow a Beta(1,10) distribution.",
    sep="\n"
  )
)

