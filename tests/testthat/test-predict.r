context("Prediction functions")

test_that("Prediction works", {
  
  # Jmorr tree
  tree  <- matrix(c(1, 2, 1, 3, 3, 4, 3, 5), ncol=2, byrow = TRUE)
  tree  <- as.phylo(tree)
  X     <- c(0, 0, 0)
  atree <- new_aphylo(tip.annotation = X, tree = tree)
  
  # Model parameters
  Pi  <- .3
  psi <- c(.01, .02)
  mu  <- c(.2, .1)
  eta <- c(1, 1)
  
  # John's spreadsheet
  ans0 <- c(
    0.009268032907,
    0.006404403815,
    0.006831065254,
    0.006041731783,
    0.006041731783
  )[with(atree$tree, c(tip.label, node.label))]
  
  # Computing brute-force and using the pre-order
  ans1 <- predict_brute_force(atree, psi, mu_d = mu, mu_s = mu, Pi)
  ans2 <- predict_pre_order(atree, psi, mu_d = mu, mu_s = mu, eta, Pi)
  
  
  expect_equivalent(ans0, ans1$posterior)
  expect_equivalent(ans1$posterior, ans2)
  
  # Random test
  set.seed(122331)
  atree <- raphylo(4, psi = psi, mu_d = mu, mu_s = mu, eta = eta, Pi = Pi)
  
  ans0 <- predict_brute_force(atree, psi, mu_d = mu, mu_s = mu, Pi)
  ans1 <- predict_pre_order(atree, psi, mu_d = mu, mu_s = mu, eta, Pi)
  
  expect_equivalent(ans0$posterior, ans1[,1])
  
})

test_that("Calling the prediction function works", {
  
  set.seed(137245)
  
  x <- raphylo(10)
  x_obs <- rdrop_annotations(x, .5)
  res   <- suppressWarnings(aphylo_mcmc(x_obs ~ psi + mu_d + Pi, priors = bprior()))
  
  ans0 <- predict_pre_order(
    x    = x_obs,
    psi  = res$par[c("psi0", "psi1")],
    mu_d = res$par[c("mu_d0", "mu_d1")],
    mu_s = res$par[c("mu_d0", "mu_d1")],
    eta  = c(1,1)/2, #res$par[c("eta0", "eta1")],
    Pi   = res$par["Pi"]
    )
  
  ans1 <- predict(res)
  
  expect_silent(plot(prediction_score(res)))
  expect_equivalent(ans0, ans1)
  
})

# Prediction score function ----------------------------------------------------

test_that("Best vs Worse Prediction score", {
  
  set.seed(123)
  y <- sample(c(0,1), 20, replace = TRUE)
  
  # Perfect prediction score
  ans0 <- prediction_score(cbind(y), cbind(y))
  expect_equivalent(ans0$obs/ans0$worse, 0)
  
  # Worse
  ans1 <- prediction_score(cbind(y), 1 - cbind(y))
  expect_equivalent(ans1$obs/ans1$worse, 1)
  
})

test_that("Random prediction score", {
  
  set.seed(123)
  y <- cbind(sample(c(0,1), 20, replace = TRUE, prob = c(.8, .2)))
  a <- .3
  
  p0 <- mean(aphylo:::predict_random(1, y, diag(20), alpha = a))
  p1 <- aphylo:::prediction_score_rand(y, diag(20), a)
  
  expect_equivalent(p0, p1, tol = 1e-1)
  
})

