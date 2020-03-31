set.seed(21354131)
x1 <- raphylo(20, P = 2)

expect_equal(x1[, 1]$tree, x1[, 2]$tree)
expect_equal(x1[, 1]$offspring, x1[, 2]$offspring)

expect_equal(
  x1[c(1, 4, 10),],
  rbind(x1$tip.annotation, x1$node.annotation)[c(1, 4, 10),]
  )

expect_equal(
  x1[c(1, 4, 10), 1],
  rbind(x1$tip.annotation, x1$node.annotation)[c(1, 4, 10), 1, drop=FALSE]
  )

expect_error(x1[], "specify")
expect_equal(
  x1[, 2, drop = TRUE],
  rbind(x1$tip.annotation[,2,drop=FALSE], x1$node.annotation[,2,drop=FALSE])
  )
expect_warning(x1[1,,drop = TRUE])

expect_error({
  x1[1,] <- c(1, 5)
}, "5")

expect_output(c(x1, x1), "2 annotated")
