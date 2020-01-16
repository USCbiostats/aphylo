set.seed(21512)
x <- raphylo(20)

x$tip.annotation[] <- rep(c(0,1), 10)
expect_equal(balance_ann(x), 1)

x$tip.annotation[] <- rep(0, 20)
expect_equal(balance_ann(x), 0)

x$tip.annotation[] <- rep(1, 20)
expect_equal(balance_ann(x), 0)

expect_error(balance_ann(x$tree))
