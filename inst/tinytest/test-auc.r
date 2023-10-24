
set.seed(1231)
x <- rnorm(100)
y <- as.integer((.2*x + rnorm(100)) > 0)
p <- stats::predict(stats::glm(y~0+x, family=binomial("probit")), type="response")


ans0 <- auc(p, y, 100)
ans1 <- AUC::auc(AUC::roc(p, as.factor(y)))

# Default way
expect_equal(ans0$auc, ans1, tol=0.01)

# Now using lists
pscore0 <- prediction_score(
  as.list(x),
  as.list(y)
)

# Now using vectors
pscore1 <- prediction_score(
  x,
  y
)

pscore2 <- # Now using matrix
prediction_score(
  cbind(x),
  cbind(y)
)

pscore3 <- # Now using a data frame
prediction_score(
  data.frame(x),
  data.frame(y)
)

expect_equal(
  pscore0$auc,
  pscore1$auc
)

expect_equal(
  pscore0$auc,
  pscore2$auc
)

expect_equal(
  pscore0$auc,
  pscore3$auc
)

