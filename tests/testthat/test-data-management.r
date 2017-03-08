context("Checking data-management methods")

# As phylo methods -------------------------------------------------------------
test_that("As phylo conversion and methods", {
  data(experiment)
  data(tree)
  
  ans <- get_offspring(
    experiment, "LeafId", 
    tree, "NodeId", "ParentId"
  )
  
  expect_s3_class(as.phylo(ans), "phylo")
  expect_is(plot(ans), "list")
  expect_output(summary(as.phylo(ans)), "Phylogenetic tree")
})