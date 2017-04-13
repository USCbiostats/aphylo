context("Checking data-management methods")

# As phylo methods -------------------------------------------------------------
test_that("As phylo conversion and methods", {
  data(fakeexperiment)
  data(faketree)
  
  ans <- new_aphylo(fakeexperiment, faketree[,c("ParentId", "NodeId")], "LeafId")
  
  expect_s3_class(as.phylo(ans), "phylo")
  expect_is(plot(ans), "list")
  expect_output(summary(as.phylo(ans)), "Phylogenetic tree")
})