# context("Area Under the Curve")

# test_that("Coincides with what the AUC::auc(roc()) function returns.", {
  
  set.seed(1231)
  x <- rnorm(100)
  y <- as.integer((.2*x + rnorm(100)) > 0)
  p <- stats::predict(stats::glm(y~0+x, family=binomial("probit")), type="response")
  
  
  ans0 <- auc(p, y, 100)
  ans1 <- AUC::auc(AUC::roc(p, as.factor(y)))
  
  expect_equal(ans0$auc, ans1, tol=0.01)
  
# })


