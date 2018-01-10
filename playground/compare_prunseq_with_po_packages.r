# This file compares the `original-w-po` version of `aphylo` with the new-using-pseq
# version of it in which. The using-po uses the pruning sequence as generated
# by ape::postorder. Moreover, it removes all the po-related classes and
# methods as these are no longer required. A simplier version of the package

rm(list = ls())

# Generating data
set.seed(1)
n <- 20
tree <- ape::rtree(n)
ann  <- data.frame(fun = sample(c(0,9,1), n, TRUE))

# New estimates
system("git checkout loglike-using-po")
devtools::load_all()

dat_new  <- new_aphylo(tree=tree, tip.annotation = ann)
ans_new <- aphylo_mle(
  dat_new,
  par     = rep(.01, 5),
  lower   = 1e-4,
  upper   = 1 - 1e-4
  )

# Old estimates
system("git checkout original-w-po")
devtools::unload("aphylo")
devtools::load_all()

dat_old  <- new_aphylo(
  edges = tree$edge,
  annotations = cbind(1:(2*n-1),rbind(ann, cbind(fun=rep(9, tree$Nnode))))
  )
ans_old <- aphylo_mle(
  dat_old,
  par     = rep(.01, 5),
  lower   = 1e-4,
  upper = 1 - 1e-4
  )
