# context("Read panther")

# test_that("Reading panther data works", {
  path <- system.file("tree.tree", package="aphylo")
  ans  <- read_panther(path)
  
  expect_equal(class(ans$tree), "phylo")
  expect_equal(nrow(ans$internal_nodes_annotations), Nnode(ans$tree))
# })