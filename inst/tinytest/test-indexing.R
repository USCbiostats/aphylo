set.seed(21354131)
x1 <- raphylo(20, P = 2)

expect_equal(x1[, 1]$tree, x1[, 2]$tree)
expect_equal(x1[, 1]$offspring, x1[, 2]$offspring)

expect_equal(x1[c(1, 4, 10),], rbind(x1$tip.annotation, x1$node.annotation)[c(1, 4, 10),])
expect_equal(x1[c(1, 4, 10), 1], rbind(x1$tip.annotation, x1$node.annotation)[c(1, 4, 10), 1, drop=FALSE])
