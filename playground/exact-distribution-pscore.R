library(aphylo)

set.seed(12444)
atree <- raphylo(20)
plot(atree)

model  <- aphylo_mcmc(atree ~ psi + mu + Pi, priors = bprior())
pscore <- prediction_score(model)

pscore_rand_dist <- aphylo:::predict_random(
  P     = 1L,
  A     = pscore$expected[1:Ntip(model)],
  G_inv = diag(length(pscore$predicted[1:Ntip(model)])),
  alpha = .8, R = 1e5L
  )

# Function to calculate the probability of obtaining k
prob_pscore <- function(x, n0, n1, alpha) {
  
  ans <- sum(vapply(0:x, function(k) {
    sum(dbinom(0:k, n1, 1 - alpha)*dbinom(k - 0:k, n0, alpha))
  }, double(1)))
  
  if (ans > .5)
    (1 - ans)
  else 
    ans
  
}

n0 <- sum(pscore$expected[1:Ntip(model)] == 0)
n1 <- sum(pscore$expected[1:Ntip(model)] == 1)
prob_pscore(pscore$random, n0, n1, .8)

mean(pscore_rand_dist <= pscore$random)
