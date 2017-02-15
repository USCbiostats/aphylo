context("MCMC")

test_that("Reasonable values", {
  
  # Simulating data
  set.seed(981)
  
  D <- rnorm(1000, 0)
  
  # Preparing function
  fun <- function(x) {
    res <- log(dnorm(D, x))
    if (any(is.infinite(res) | is.nan(res)))
      return(.Machine$double.xmin)
    sum(res)
  }
  
  # Running the algorithm and checking expectation
  ans <- suppressWarnings(
    MCMC(fun, 1, 1e4, burnin = 1000, ub = 3, lb = -3, scale = 1)
  )
  expect_equal(mean(ans), 0, tolerance = 0.1, scale = 1)
})