rm(list = ls())

library(aphylo)

n <- 200
N <- 500

set.seed(111222)
PAR <- lapply(1:N, function(i) rbeta(5, 2, 20))
DAT <- lapply(PAR, function(p) sim_annotated_tree(n, psi=p[1:2], mu=p[3:4], Pi=p[5]))

plot_LogLike(DAT[[1]], psi = c(.1, .1), mu = c(.1, .1), Pi=.1)

est <- function(d) {
  tryCatch({
    aphylo_mle(params = rep(.05,5), dat = d, #, control=list(nbatch=5e3, burnin=1e3, thin=10),
                priors = function(u) dbeta(u, 2, 10))
  }, error = function(e) e)
}

cl <- parallel::makeForkCluster(4)
ANS <- parallel::parLapply(cl, DAT, est)
# ANS <- lapply(DAT, est)
parallel::stopCluster(cl)

PARest <- lapply(ANS, "[[", "par")

OK <- which(grepl("aphy", sapply(sapply(ANS, class), "[[", 1)))

PARest <- do.call(rbind, PARest[OK])
PARh0  <- do.call(rbind, lapply(1:length(OK), function(i) rbeta(5, 2, 20)))

bias   <- PARest - do.call(rbind, PAR[OK])
biash0 <- PARh0 - do.call(rbind, PAR[OK])

oldpar <- par(no.readonly = TRUE)
par(mfrow = c(1, 2))
boxplot(bias, main = "Non random", ylim = c(-.5, .5))
boxplot(biash0, main = "Random", ylim = c(-.5, .5))
par(oldpar)
summary(bias)
summary(biash0)

hist(sapply(ANS[OK], "[[", "counts"))
prop.table(table(abs(bias) < abs(biash0)))
