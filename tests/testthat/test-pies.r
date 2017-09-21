context("Pies work")

test_that("We can only pie numbers", {
  expect_error(piechart(letters), "must be numeric")
  expect_silent(piechart(1:10, labels = letters[1:10]))
  expect_silent(piechart(1:10, labels = letters[1:10], radius = 1:10 /10,
                         slice.off = 1:10 /20, col=colors()[1:10]))
})