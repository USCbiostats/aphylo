rm(list = ls())

tree <- read.table("data-raw/pthr11848.dag.txt",
                  col.names = c("NodeId", "TypeId", "ParentId"))
tree <- as.matrix(tree)
save(tree, file = "data/tree.rda")
