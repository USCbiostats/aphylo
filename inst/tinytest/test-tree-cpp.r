# context("Testing tree.cpp")

# test_that("Fast table works", {
  
  set.seed(1)
  dat <- cbind(sample.int(10, 500, TRUE))
  ans0 <- table(dat)
  ans1 <- aphylo:::fast_table(dat)
  ans2 <- aphylo:::fast_table_using_labels(dat, 1L:10L)
  
  expect_equivalent(as.vector(ans0), ans1[,2])
  expect_equivalent(as.vector(ans0), ans2)
  
  ans3 <- aphylo:::fast_table_using_labels(dat, c(1, 4, 10))
  expect_equivalent(as.vector(ans0)[c(1, 4, 10)], ans3)
# })