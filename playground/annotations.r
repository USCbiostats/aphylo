rm(list = ls())

library(aphylo)
library(ape)

# Reading in a TREE and annotations
dat <- read.table("../aphylo_paper/data-raw/experimental_function", sep="\t",
                  stringsAsFactors = FALSE,
                  col.names = c("node_id", "lead_id", "go_id", "annotation")
                  )

# Getting the panther book
dat$book_id <- gsub("[:].+$", "", dat$node_id)

# Coding annotations
dat$status <- with(
  dat, 
  ifelse(annotation == "NOT", 0L, ifelse(annotation == "", NA, 1L))
  )

library(dplyr)
counts <- group_by(dat, book_id, go_id) %>%
  summarize(
    nas   = sum(is.na(status)),
    zeros = sum(status[!is.na(status)] == 0),
    ones  = sum(status[!is.na(status)] == 1)
  )

counts <- counts[with(counts, order(zeros, ones, decreasing = TRUE)),]

# This is the set of books that we should use
counts_relevant <- as.data.frame(subset(counts, (zeros != 0) & (ones != 0)))

# Book PTHR23255 and go GO:0050431 is nice
# PTHR10788
# GO:0004805
tree <- read.tree("../PANTHER11.1/books/PTHR10788/tree.tree", nmax=1)
tree$tip.label <- paste("PTHR10788", tree$tip.label, sep=":")
ids <- tree$tip.label

A <- subset(dat, (book_id == "PTHR10788") & (node_id %in% ids) & (go_id == "GO:0004805"))

# Creating the aphylo object
atree <- new_aphylo(
  data.frame(A$node_id, `GO:0004805`= A$status, stringsAsFactors = FALSE,
             check.names = FALSE),
  as_po_tree(tree)
)

table(atree$annotations)

ans_mle  <- aphylo_mle(atree)
ans_mcmc <- aphylo_mcmc(ans_mle$par, atree, priors = function(x) dbeta(x, 2, 30),
                        control=list(nbatch=1e5, thin=100, burnin=1e4, nchains=5))


cbind(
  predict(ans_mcmc, what = names(which(ans_mcmc$dat$annotations[,1] != 9))),
  ans_mcmc$dat$annotations[names(which(ans_mcmc$dat$annotations[,1] != 9)),]
  )

cbind(
  predict(ans_mle, what = names(which(ans_mle$dat$annotations[,1] != 9))),
  ans_mle$dat$annotations[names(which(ans_mle$dat$annotations[,1] != 9)),]
)
