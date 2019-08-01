# context("Geodesic distances")

# test_that("geodesic", {
  # Simple (obvious) test
  edges <- matrix(c(0,1,1,2,2,3,3,4), ncol=2, byrow = TRUE)
  ans0 <- matrix(c(0:4, 1, 0:3, 2:0, 1,2, 3:0, 1, 4:0), ncol=5)
  ans1 <- approx_geodesic(edges + 1)
  expect_equal(ans0, ans1)
  
  # A more complicated test
  edges <- matrix(c(3,6,3,8,2,7,2,4,1,5,1,2,0,3,0,1), byrow = TRUE, ncol=2)
  ans0 <- matrix(0, ncol=9, nrow = 9)
  ans0[1, ]         <- c(0, 1, 2, 1, 3, 2, 2, 3, 2)
  ans0[2, -1]       <- c(0, 1, 2, 2, 1, 3, 2, 3)
  ans0[3, -c(1, 2)] <- c(0, 3, 1, 2, 4, 1, 4)
  ans0[4, -c(1:3)]  <- c(0, 4, 3, 1, 4, 1)
  ans0[5, -c(1:4)]  <- c(0, 3, 5, 2, 5)
  ans0[6, -c(1:5)]  <- c(0, 4, 3, 4)
  ans0[7, -c(1:6)]  <- c(0, 5, 2)
  ans0[8, -c(1:7)]  <- c(0, 5)
  ans0 <- ans0 + t(ans0)
  ans1 <- approx_geodesic(edges + 1)
  
  expect_equal(ans0, ans1)
# })