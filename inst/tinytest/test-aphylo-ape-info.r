set.seed(18231)

x <- c(raphylo(20), raphylo(10, P = 2))
y <- c(x[[1]]$tree, x[[2]]$tree)
z <- new_aphylo_pruner(x)


res0 <- aphylo_mle(x ~ psi + mu_d + mu_s + Pi)

expect_equal(Ntip(x), c(20, 10))
expect_equal(Ntip(y), c(20, 10))
expect_equal(Ntip(res0), c(20, 10))
expect_equal(Ntip(z), c(20, 10))

expect_equal(Nnode(x), c(19, 9))
expect_equal(Nnode(y), c(19, 9))
expect_equal(Nnode(res0), c(19, 9))
expect_equal(Nnode(z), c(19, 9))

expect_equal(Nann(x), c(1, 2))
expect_equal(Nann(res0), c(1, 2))
expect_equal(Nann(z), c(1, 2))

expect_equal(Nannotated(x), c(20, 10))
expect_equal(Nannotated(res0), c(20, 10))
expect_equal(Nannotated(z), c(20, 10))

expect_equal(Ntrees(x), 2)
expect_equal(Ntrees(y), 2)
expect_equal(Ntrees(res0), 2)
expect_equal(Ntrees(z), 2)


