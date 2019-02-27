context("Prediction functions")

test_that("Prediction works", {
  
  # Jmorr tree
  tree  <- matrix(c(1, 2, 1, 3, 3, 4, 3, 5), ncol=2, byrow = TRUE)
  tree  <- as.phylo(tree)
  X     <- c(0, 0, 0)
  atree <- new_aphylo(X, tree)
  
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
  ans1 <- predict_brute_force(atree, psi, mu, Pi)
  ans2 <- predict_pre_order(atree, psi, mu, eta, Pi)
  
  
  expect_equivalent(ans0, ans1$posterior)
  expect_equivalent(ans1$posterior, ans2)
  
  # Random test
  set.seed(122331)
  atree <- sim_annotated_tree(4, psi = psi, mu = mu, eta = eta, Pi = Pi)
  
  ans0 <- predict_brute_force(atree, psi, mu, Pi)
  ans1 <- predict_pre_order(atree, psi, mu, eta, Pi)
  
  expect_equivalent(ans0$posterior, ans1[,1])
  
})

test_that("Calling the prediction function works", {
  
  set.seed(137245)
  
  x <- sim_annotated_tree(10)
  x_obs <- rdrop_annotations(x, .5)
  res   <- suppressWarnings(aphylo_mcmc(x_obs ~ psi + mu + Pi, priors = bprior()))
  
  ans0 <- predict_pre_order(
    x   = x_obs,
    psi = res$par[c("psi0", "psi1")],
    mu  = res$par[c("mu0", "mu1")],
    eta = c(1,1)/2, #res$par[c("eta0", "eta1")],
    Pi  = res$par["Pi"]
    )
  
  ans1 <- predict(res)
  
  expect_silent(plot(prediction_score(res)))
  expect_equivalent(ans0, ans1)
  
})
