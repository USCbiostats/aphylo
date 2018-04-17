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
  expect_equivalent(ans1$posterior, ans1$posterior)
  
  # Random test
  set.seed(122331)
  atree <- sim_annotated_tree(6, psi = psi, mu = mu, eta = eta, Pi = Pi)
  
  ans0 <- predict_brute_force(atree, psi, mu, Pi)
  ans1 <- predict_pre_order(atree, psi, mu, eta, Pi)
  
  all.equal(ans0$posterior, ans1$posterior[,1])
  
})