rm(list = ls())

library(mcmc)

# Parameters
set.seed(1231)
n <- 1e3
pars <- c(mean = 2.6, sd = 3)

# Generating data and writing the log likelihood function
D <- rnorm(n, pars[1], pars[2])
fun <- function(x) sum(log(dnorm(D, x[1], x[2])))

# Running the MCMC
ans <- metrop(fun, c(.5,.5), nbatch = 1e5, scale = .005)

# Removing the first 5000
ans$batch <- ans$batch[-c(1:5000),]

# Checking out the results
graphics.off()
png("playground/testing_metrop_in_mcmc.png", height = 800, width = 600)
oldpar <- par(no.readonly = TRUE)
par(mfrow = c(2,1))
boxplot(ans$batch, 
        main = expression("Posterior distribution of"~mu~and~sigma),
        names =  expression(mu, sigma), horizontal = TRUE,
        col  = blues9[c(4,9)],
        sub = bquote(mu == .(pars[1])~", and"~sigma == .(pars[2]))
        )
abline(v = par, col  = blues9[c(4,9)], lwd = 2, lty = 2)

plot(apply(ans$batch, 1, fun), type = "l",
     main = "LogLikelihood",
     ylab = expression(L("{"~mu,sigma~"}"~"|"~D)) 
     )
par(oldpar)
dev.off()