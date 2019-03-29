library(aphylo)

rm(list = ls())

# Checking whether reducing the peeling sequence helps
set.seed(1)
x <- raphylo(400)
x$node.annotation[] <- 9
x$tip.annotation[sample.int(400, 300)] <- 9
x0 <- x
x1 <- x
x1$pseq <- reduce_pseq(x$pseq, with(x, rbind(tip.annotation, node.annotation)), x$offspring)

microbenchmark::microbenchmark(
  LogLike(x0, c(.1,.1), c(.2,.2), .1),
  LogLike(x1, c(.1,.1), c(.2,.2), .1), times = 5e4,
  unit = "relative"
)


LogLike(x0, c(.1,.1), c(.2,.2), .1)$ll
LogLike(x1, c(.1,.1), c(.2,.2), .1)$ll

length(x0$pseq)
length(x1$pseq)

length(x1$pseq)/length(x0$pseq)
