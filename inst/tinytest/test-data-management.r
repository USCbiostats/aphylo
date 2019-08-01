# context("Checking data-management methods")

data(fakeexperiment)
data(faketree)

# As phylo methods -------------------------------------------------------------
# test_that("As phylo conversion and methods", {
  
  tree <- as.phylo(faketree)
  ans <- new_aphylo(tip.annotation = fakeexperiment[,-1], tree = tree)
  
  expect_silent(print(ans))
  expect_equal(class(as.phylo(ans)), "phylo")
  expect_true(is.null(plot(ans)))
  expect_silent(summary(as.phylo(ans)))
  expect_silent(summary(raphylo(10, P=4)))
# })

# Conversion -------------------------------------------------------------------
# test_that("Can return to the original labeling", {
  ans0  <- faketree
  phylo <- as.phylo(ans0)
  ans1  <- phylo$edge
  
  
  ans1[] <- with(phylo, c(tip.label, node.label))[ans1[]]
  
  expect_equivalent(ans0, ans1)
# })


