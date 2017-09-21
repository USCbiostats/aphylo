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
  
  expect_equivalent(ds, degseq(sim_tree(n)))
})