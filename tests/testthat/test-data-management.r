context("Checking data-management methods")

data(fakeexperiment)
data(faketree)

# As phylo methods -------------------------------------------------------------
test_that("As phylo conversion and methods", {
  ans <- new_aphylo(fakeexperiment, faketree[,c("ParentId", "NodeId")], "LeafId")
  
  expect_s3_class(as.apephylo(ans), "phylo")
  expect_is(plot(ans), "list")
  expect_output(summary(as.apephylo(ans)), "Phylogenetic tree")
})

# Conversion -------------------------------------------------------------------
test_that("Can return to the original labeling", {
  ans0 <- faketree
  ans1 <- as_ape_tree(ans0)
  
  ans1[] <- attr(ans1, "labels")[ans1[]]
  
  expect_equivalent(ans0, ans1)
})