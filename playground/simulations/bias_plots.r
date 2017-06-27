rm(list = ls())

library(ggplot2)
library(magrittr)

load("playground/simulations/bias.rda")

# Reshaping
bias_vars <- paste0(c("psi0", "psi1", "mu0", "mu1", "Pi"), "_bias")
info_vars <- c("index", "Method", "size_tag", "miss_tag", "PropLeafs_tag")
bias <- do.call(rbind, lapply(bias_vars, function(x) {
  ans <- data.frame(parameter = x, bias[,c(info_vars, x),drop=FALSE])
  colnames(ans)[7] <- "bias"
  ans
}
))

#+ plotting, echo=TRUE
# Plot -------------------------------------------------------------------------
graphics.off()
sizelvls <- levels(bias$size_tag)
for (i in 1:4) {
  # Clearing plot space and creating the pdf
  
  pdf(sprintf("playground/simulations/bias_trees_of_size_%s.pdf", sizelvls[i]))
  
  # Nobservations in this group
  nobs <- bias %>% dplyr::filter(as.integer(size_tag) == i) %>%
    subset(select=index) %>% unique %>% nrow
  
  
  # Creating the plot
  p <- bias %>% dplyr::filter(as.integer(size_tag) == i) %>%
    # Creating the boxplot
    ggplot(aes(parameter, bias)) + geom_boxplot(aes(colour = Method)) +
    
    # Adding an horizontal line, and spliting my % of missings
    geom_hline(yintercept = 0, lty=2) + facet_grid(miss_tag ~ .) +
    
    # Adding a title
    ggtitle(
      label    = sprintf("Bias distribution for trees of size %s", levels(bias$size_tag)[i]),
      subtitle = sprintf("# of observations: %i", nobs)
    ) +
    
    ylim(-.75,.75)
  
  # Printing it on screen (need to do that explicitly on a loop)
  print(p)
  message("Level ", levels(bias$size_tag)[i], " done.")
  
  dev.off()
  
}

# Computing table of biases
bias_mat <- with(bias, array(
  dim = c(nlevels(miss_tag), nlevels(Method), nlevels(parameter)),
  dimnames = list(
    levels(miss_tag),
    levels(Method),
    levels(parameter)
  ))
)

for (m in levels(bias$miss_tag)) 
  for (me in levels(bias$Method)) 
    for (p in levels(bias$parameter))
      bias_mat[m, me, p] <- mean(subset(bias, miss_tag == m & Method == me & parameter == p)$bias)


var_mat <- with(bias, array(
  dim = c(nlevels(miss_tag), nlevels(Method), nlevels(parameter)),
  dimnames = list(
    levels(miss_tag),
    levels(Method),
    levels(parameter)
  ))
)

for (m in levels(bias$miss_tag)) 
  for (me in levels(bias$Method)) 
    for (p in levels(bias$parameter))
      var_mat[m, me, p] <- var(subset(bias, miss_tag == m & Method == me & parameter == p)$bias)

