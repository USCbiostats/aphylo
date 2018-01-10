context("Simulation of Annotated Phylogenetic Trees")

# Function to compute degree distribution
degseq <- function(x) {
  ind <- as.data.frame(table(x[,2]), responseName = "ind")
  oud <- as.data.frame(table(x[,1]), responseName = "oud")
  
  ans <- merge(ind, oud, by="Var1", all = TRUE)
  ans[is.na(ans)] <- 0
  as.vector(table(ans$ind + ans$oud))
}

test_that("Simulating Trees", {
  set.seed(121)
  n <- 100
  d <- vector("integer", n)
  
  ds <- c(n, 1, n-2)
  
  expect_equivalent(ds, degseq(sim_tree(n)$edge))
})

checkout_annotations <- function(x) {
  sapply(
    lapply(x, "[[", "tip.annotation"),
    aphylo:::fast_table_using_labels, ids=c(0,1)
  )
}

test_that("Simulating informative annotated trees", {
  
  set.seed(1)
  ans <- lapply(1:100, function(i) sim_annotated_tree(5))
  ans <- checkout_annotations(ans)
  
  # There must be zeros and ones in all trees!
  expect_true(all(ans[1,] >= 1))
  expect_true(all(ans[2,] >= 1))
  
})

test_that("Dropping annotations", {
  set.seed(1)
  sizes <- ceiling(runif(1e3, .01)*200)
  ans <- lapply(sizes, sim_annotated_tree, Pi = .5, psi=c(0,0), mu=c(.5, .5))
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
  
  expect_equal(mean(ans3), 1/3, tol = .1)
  
})