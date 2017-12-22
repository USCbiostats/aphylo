context("Checking data-management methods")

data(fakeexperiment)
data(faketree)

# As phylo methods -------------------------------------------------------------
test_that("As phylo conversion and methods", {
  ans <- new_aphylo(tip.annotation = fakeexperiment[,-1], tree = faketree)
  
  expect_s3_class(as.phylo(ans), "phylo")
  expect_s3_class(plot(ans), "ggplot")
  expect_output(summary(as.phylo(ans)), "Phylogenetic tree")
  expect_output(summary(sim_annotated_tree(10, P=4)), "Distri")
})

# Conversion -------------------------------------------------------------------
test_that("Can return to the original labeling", {
  ans0  <- faketree
  phylo <- as.phylo(ans0)
  ans1  <- phylo$edge
  
  
  ans1[] <- with(phylo, c(tip.label, node.label))[ans1[]]
  
  expect_equivalent(ans0, ans1)
})


test_that("Listing leafs", {
  set.seed(1)
  dat0 <- sim_tree(100)
  dat1 <- as.phylo(dat0)
  
  dat0 <- leafs(dat0)
  dat1 <- leafs(dat1)
  
  expect_true(all(sort(dat0) == sort(dat1)))
})