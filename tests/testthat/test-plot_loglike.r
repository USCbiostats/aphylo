context("plot_LogLike")

test_that("plotll", {
  set.seed(1)
  dat <- sim_annotated_tree(20)
  
  # aphylo method
  expect_silent(plot_LogLike(dat))

  # phylo_mle method
  ans <- phylo_mle(dat)
  expect_silent(plot_LogLike(ans))
  
})