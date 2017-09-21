context("Checking data-management methods")

data(fakeexperiment)
data(faketree)

# As phylo methods -------------------------------------------------------------
test_that("As phylo conversion and methods", {
  ans <- new_aphylo(fakeexperiment, faketree)
  
  expect_s3_class(as.apephylo(ans), "phylo")
  expect_s3_class(plot(ans), "ggplot")
  expect_output(summary(as.apephylo(ans)), "Phylogenetic tree")
})

# Conversion -------------------------------------------------------------------
test_that("Can return to the original labeling", {
  ans0 <- faketree
  ans1 <- as_ape_tree(ans0)
  
  ans1[] <- attr(ans1, "labels")[ans1[]]
  
  expect_equivalent(ans0, ans1)
})

# ------------------------------------------------------------------------------
test_that("as_po tree", {
  
  # Generating a random tree
  set.seed(100)
  ans0 <- sim_tree(50)$edges
  
  attr(ans0, "labels") <-
    structure(
      paste0(1:1000,letters)[1:length(attr(ans0, "labels"))],
      names = names(attr(ans0, "labels")))
  
  ans1 <- as_po_tree(as.apephylo(ans0))
  
  # Plotting should work fine
  expect_silent(plot(ans0, show.node.label=TRUE))
  
  # Replacing labels
  ans0[] <- attr(ans0, "labels")[ans0[]+1]
  ans1[] <- attr(ans1, "labels")[ans1[]+1]
  
  ans0 <- ans0[order(ans0[,1], ans0[,2]), ]
  ans1 <- ans1[order(ans1[,1], ans1[,2]), ]
  
  # These should be the same afterwards
  expect_equal(ans0, ans1)
  
  
})

test_that("Listing leafs", {
  set.seed(1)
  dat0 <- sim_tree(100)
  dat1 <- as.apephylo(dat0$edges)
  
  dat0 <- leafs(dat0$edges)
  dat1 <- leafs(dat1)
  
  expect_true(all(sort(dat0) == sort(dat1)))
})