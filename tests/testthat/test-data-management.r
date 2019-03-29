context("Checking data-management methods")

data(fakeexperiment)
data(faketree)

# As phylo methods -------------------------------------------------------------
test_that("As phylo conversion and methods", {
  ans <- new_aphylo(tip.annotation = fakeexperiment[,-1], tree = faketree)
  
  expect_output(print(ans), "\\(leafs\\) annotations")
  expect_s3_class(as.phylo(ans), "phylo")
  expect_null(plot(ans))
  expect_output(summary(as.phylo(ans)), "Phylogenetic tree")
  expect_output(summary(raphylo(10, P=4)), "Distri")
})

# Conversion -------------------------------------------------------------------
test_that("Can return to the original labeling", {
  ans0  <- faketree
  phylo <- as.phylo(ans0)
  ans1  <- phylo$edge
  
  
  ans1[] <- with(phylo, c(tip.label, node.label))[ans1[]]
  
  expect_equivalent(ans0, ans1)
})


