set.seed(1231512)

x <- raphylo(20)
ans0 <- aphylo_cv(
  x ~ psi + mu_d + mu_s + Pi, 
  control = list(nsteps = 500, nchains = 1, burnin = 0)
  )

ans1 <- aphylo_mcmc(
  x ~ psi + mu_d + mu_s + Pi, 
  control = list(nsteps = 500, nchains = 1, burnin = 0)
  )


expect_equal(coef(ans0$estimates), coef(ans1), tol = .2)
expect_output(print(ans0$auc), "AUC")
expect_silent(plot(ans0$auc))

x <- replicate(2, raphylo(20), simplify = FALSE)
x <- do.call(c, x)
ans0 <- aphylo_cv(
  x ~ psi + mu_d + mu_s + Pi, 
  control = list(nsteps = 500, nchains = 1, burnin = 0)
)

ans1 <- aphylo_mcmc(
  x ~ psi + mu_d + mu_s + Pi, 
  control = list(nsteps = 500, nchains = 1, burnin = 0)
)
