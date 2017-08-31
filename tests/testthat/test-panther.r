context("Read panther")

test_that("Reading panther data works", {
  path <- system.file("tree.tree", package="aphylo")
  ans  <- read_panther(path)
  
  expect_output(print(ans), "Phylogenetic tree with 145 tips and 107 internal nodes.")
})