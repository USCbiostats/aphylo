# context("plot_LogLike")

# test_that("plotll", {
  set.seed(1)
  dat <- raphylo(20)
  
  # aphylo method
  expect_true(is.null(plot_logLik(dat)))

  # phylo_mle method
  ans <- suppressWarnings(aphylo_mle(dat ~ mu_d + psi))
  expect_true(is.null(plot_logLik(ans)))
  
# })
