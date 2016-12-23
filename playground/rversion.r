rm(list=ls())

# Reading data
dag <- read.table("data/pthr11848.dag.txt", sep="\t", header = FALSE,
                  col.names = c("NodeId", "TypeId", "ParentId"))
dat <- read.table("data/pthr11848_sorted.txt", sep = "\t", header = FALSE)
colnames(dat) <- c(sprintf("f%02d",1:(ncol(dat)-1)), "LeafId")

# Is it already sorted?
stopifnot(
  all(dag$NodeId == sort(dag$NodeId)) & all(dat$LeafId == sort(dat$LeafId))
  )

# Checking for which there is data
dat_notin_dag <- which(!(dat$LeafId %in% c(dag$NodeId, dag$ParentId)))
dag_notin_dat <- which(!(c(dag$NodeId, dag$ParentId) %in% dat$LeafId))

# Filling dat with NAs if NodeId is not in LeafId
dat[dag_notin_dat,1:3] <- NA

length(dat_notin_dag)
length(dag_notin_dat)

# 
N <- nrow(dag)
parent <- vector("integer",N)
parent[] <- NA
parent[1] <- -1L
for (n in 2:N) {
  found <- 0L
  par   <- 0L
  while (!found) {
    if (with(dag, ParentId[n] == NodeId[par+1L])) {
      parent[n] <- par
      found     <- 1L
    }
    par <- par + 1L
    if ((par + 1L) > N) break
  }
}

OffId <- lapply(dag$NodeId, function(x) {
  which(dag$ParentId == x)
})
