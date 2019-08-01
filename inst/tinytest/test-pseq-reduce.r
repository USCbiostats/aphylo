# context("Reducing pruning sequence")

# test_that("Same results", {
  
  set.seed(1235541)
  x <- raphylo(100)
  x <- rdrop_annotations(x, .5)
  
  control. <- list(multicore = FALSE, nchains=1, nsteps=1000, burnin=0)
  
  set.seed(1)
  ans0 <- aphylo_mcmc(x ~ psi + Pi + mu_d, priors = bprior(), reduced_pseq = FALSE,
                      control = control.)
  set.seed(1)
  ans1 <- aphylo_mcmc(x ~ psi + Pi + mu_d, priors = bprior(), reduced_pseq = TRUE,
                      control = control.)
  
  expect_identical(ans0$hist, ans1$hist)
  
  set.seed(1)
  ans0 <- aphylo_mcmc(x ~ psi + Pi + mu_d + eta, priors = bprior(), reduced_pseq = FALSE,
                      control = control.)
  set.seed(1)
  ans1 <- aphylo_mcmc(x ~ psi + Pi + mu_d + eta, priors = bprior(), reduced_pseq = TRUE,
                      control = control.)
  
  expect_identical(ans0$hist, ans1$hist)
  
# })

# pseq1 <- x$pseq
# pseq2 <- aphylo:::reduce_pseq(
#   pseq1, 
#   with(x, rbind(tip.annotation, node.annotation)),
#   x$offspring
#   )
# 
# 
# all(pseq2 %in% pseq1)

