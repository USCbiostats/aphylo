# context("Simulation of Annotated Phylogenetic Trees")

# Function to compute degree distribution
degseq <- function(x) {
  ind <- as.data.frame(table(x[,2]), responseName = "ind")
  oud <- as.data.frame(table(x[,1]), responseName = "oud")
  
  ans <- merge(ind, oud, by="Var1", all = TRUE)
  ans[is.na(ans)] <- 0
  as.vector(table(ans$ind + ans$oud))
}

# test_that("Simulating Trees", {
  set.seed(121)
  n <- 100
  d <- vector("integer", n)
  
  ds <- c(n, 1, n-2)
  
  expect_equivalent(ds, degseq(sim_tree(n)$edge))
# })

checkout_annotations <- function(x) {
  sapply(
    lapply(x, "[[", "tip.annotation"),
    function(z) tabulate(z+1, 2)
  )
}

# test_that("Simulating informative annotated trees", {
  
  set.seed(1)
  ans <- lapply(1:100, function(i) raphylo(5, informative = TRUE))
  ans <- checkout_annotations(ans)
  
  # There must be zeros and ones in all trees!
  expect_true(all(ans[1,] >= 1))
  expect_true(all(ans[2,] >= 1))
  
# })

# test_that("Dropping annotations", {
  set.seed(1)
  sizes <- ceiling(runif(1e3, .01)*300)
  ans <- lapply(sizes, raphylo, Pi = .5, psi=c(0,0), mu_d=c(.5, .5), mu_s=c(.5, .5))
  ans0 <- checkout_annotations(ans)
  
  # Dropping half of it and keeping informative
  ans1 <- lapply(ans, rdrop_annotations, pcent=.5, informative = TRUE)
  ans1 <- checkout_annotations(ans1)
  
  expect_equal(mean(colSums(ans1)/colSums(ans0)), .5, tol = .025)
  
  # Dropping 2/3 of it and keeping informative
  ans2 <- lapply(ans, rdrop_annotations, pcent=2/3, informative = TRUE)
  ans2 <- checkout_annotations(ans2)
  
  expect_equal(mean(colSums(ans2)/colSums(ans0)), 1/3, tol = .05)
  
  # Zeros are more likely to be droped
  ans3 <- lapply(ans, rdrop_annotations, pcent = .5, informative = TRUE,
                 prob.drop.0 = 2/3)
  ans3 <- checkout_annotations(ans3)
  ans3 <- ans3[1,]/colSums(ans3)
  
  expect_equal(abs(mean(ans3) - .3333), 0, tol = .1)
  
# })

# test_that("Mislabeling", {
  
  # Ramdom tree (fully annotated) ---------------------------------------------
  set.seed(1)
  x <- raphylo(20, P=4, Pi=.02)
  
  # All are ones/zeros
  all_ones <- mislabel(x, psi = c(1, 0))
  all_zero <- mislabel(x, psi = c(0, 1))
  
  expect_true(all(all_ones$tip.annotation == 1))
  expect_true(all(all_zero$tip.annotation == 0))
  
  # All oposite
  all_flipped <- mislabel(x, psi = c(1, 1))
  expect_equivalent(x$tip.annotation, 1 - all_flipped$tip.annotation)
  
  # Same tests but setting annotations equal to 9 ------------------------------
  
  # Dropping annotations
  x <- rdrop_annotations(x, .5)
  
  # All are ones/zeros
  all_ones <- mislabel(x, psi = c(1, 0))
  all_zero <- mislabel(x, psi = c(0, 1))

  # Number of 9s is preserved
  expect_equal(sum(all_ones$tip.annotation == 9), sum(x$tip.annotation == 9))
  expect_equal(sum(all_zero$tip.annotation == 9), sum(x$tip.annotation == 9))
    
  expect_true(all(all_ones$tip.annotation[all_ones$tip.annotation != 9] == 1))
  expect_true(all(all_zero$tip.annotation[all_zero$tip.annotation != 9] == 0))
  
  # All oposite
  all_flipped <- mislabel(x, psi = c(1, 1))
  expect_equivalent(
    x$tip.annotation[x$tip.annotation != 9],
    1 - all_flipped$tip.annotation[all_flipped$tip.annotation != 9])
  
  
  
# })

# test_that("Deterministic results including node types", {
  
  # Fake tree assuring evolution
  dat2   <- as.phylo(matrix(c(1, 2, 1, 3, 2, 4, 2, 5), ncol=2, byrow = TRUE))
  
  getann <- function(x) {
    ans <- structure(
      c(x$tip.annotation, x$node.annotation),
      names = c(x$tree$tip.label, x$tree$node.label))
    ans[order(as.integer(names(ans)))]
  }
  
  tree2  <- new_aphylo(tip.annotation=rbind(0,0,0), tree = dat2)
  ans0 <- raphylo(tree = tree2, mu_d = c(1, 0), Pi = 1, psi = c(0,0),
                  eta = c(1,1), node.type = NULL)
  ans1 <- raphylo(tree = tree2, mu_d = c(0, 1), Pi = 1, psi = c(0,0),
                  eta = c(1,1), node.type = NULL)
  ans2 <- raphylo(tree = tree2, mu_d = c(1, 1), Pi = 1, psi = c(0,0),
                  eta = c(1,1), node.type = NULL)
  
  expect_equivalent(getann(ans0), rep(1, 5))        # Only gains
  expect_equivalent(getann(ans1), c(1, 0, 0, 0, 0)) # All losses
  expect_equivalent(getann(ans2), c(1, 0, 0, 1, 1)) # Mixed
  
  types <- c(1, 0, 1, 1, 1)
  tree3 <- new_aphylo(
    tree = dat2, tip.annotation = rbind(0,0,0), 
    node.type = types[dat2$node.label],
    tip.type  = types[dat2$tip.label],
    )
  
  (ans3  <- raphylo(
    tree = tree3,
    mu_s = c(0, 0), # Nothing happens in speciation nodes (1)
    mu_d = c(1, 0), # Duplication nodes gain a function
    Pi   = 0,       # Start with no function
    eta  = c(1, 1), # Everything is reported
    psi  = c(0,0),  # Nothing is misclassified
    node.type = NULL
    ))
  
  expect_equivalent(getann(ans3), c(0, 0, 0, 1, 1))
  
  types <- c(1, 0, 1, 1, 1)
  tree3 <- new_aphylo(
    tree           = dat2,
    tip.annotation = rbind(0,0,0), 
    node.type = types[dat2$node.label],
    tip.type  = types[dat2$tip.label],
    )
  
  (ans3  <- raphylo(
    tree = tree3,
    mu_s = c(0, 0), # Nothing happens in speciation nodes (1)
    mu_d = c(0, 1), # Duplication nodes loss a function
    Pi   = 1,       # Start with a function
    eta  = c(1, 1), # Everything is reported
    psi  = c(0, 0), # Nothing is misclassified
    node.type = tree3$node.type
  ))
  
  expect_equivalent(getann(ans3), c(1, 1, 1, 0, 0))
# })

