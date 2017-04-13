rm(list = ls())

faketree <- read.csv2("data-raw/faketree.csv")[,2:1]
faketree <- as.matrix(faketree)

save(faketree, file = "data/faketree.rda")
