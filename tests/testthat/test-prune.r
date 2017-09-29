context("Tree prunning")


test_that("Basic set of tests", {
  
  # A known tree
  x <- structure(
    c(2L, 2L, 1L, 1L, 0L, 0L, 3L, 4L, 6L, 2L, 1L, 5L),
    .Dim = c(6L, 2L),
    class = c("po_tree", "matrix"),
    labels = c("0", "1", "2", "3", "4", "5", "6"),
    offspring = list(
      structure(c(1, 5), .Dim = 1:2), 
      structure(c(6, 2), .Dim = 1:2),
      structure(c(3, 4), .Dim = 1:2), 
      NULL, NULL, NULL, NULL)
  )
  
  # Error messages
  expect_error(prune(x, 0), "Root")
  expect_error(prune(x, 7), "out of range")
  expect_error(prune(x, c(1,5)), "all but the root")
  
  # Removing 2 should take out 4 and 3 (so 3 nodes)
  # This shoulod cause the ids to go from 0 to 4
  z <- prune(x, 2)
  expect_true(all(aphylo:::getlabels(z) == as.character(c(0,1,5,6))))
  expect_true(length(attr(z, "labels")) == 4)
  
  # Prunning using labels instead
  attr(x, "labels") <- letters[1:7]
  expect_equal(prune(x, 2), prune(x, "c"))
  
})