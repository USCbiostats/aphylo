#' ---
#' title: "Convergence"
#' author: "George G Vega Yon"
#' date: "`r paste('This version:', Sys.time())`"
#' output: pdf_document
#' ---
#' 

#+ setup, echo=FALSE
knitr::opts_chunk$set(echo = FALSE)

#+ data-loading, cache=TRUE
load("playground/simulations/simdata3.rda")


library(ggplot2)
library(magrittr)

# Function to measure bias
bias_calc <- function(x, par0, tree, index) {
  
  # Checking if it was able to solve it or not
  if (inherits(x, "error")) 
    return(rep(NA, 9))
  
  # Number of offspring and internal nodes
  treesize <- length(tree$noffspring)
  nleafs   <- sum(tree$noffspring > 0)
  
  structure(
    c(
      index,
      treesize,
      nleafs,
      par0[6],
      x$par - par0[1:5]
      ),
    names = c(
      "index",
      "TreeSize",
      "NLeafs",
      "Missings",
      names(x$par))
    )
  
}

# Creates nice interval tags in the form of [a,b)...[a,z] (last one closed).
# All numbers must be within the interval
interval_tags <- function(x, marks) {
  
  # Generating labels
  n <- length(marks)
  l <- c(sprintf("[%.1f, %.1f)", marks[-n][-(n-1)], marks[-n][-1]),
    sprintf("[%.1f, %.1f]", marks[n-1], marks[n])
  )
  
  # Finding intervals
  x <- findInterval(x, marks, rightmost.closed = TRUE)
  factor(x, levels = 1:length(l), labels = l)
  
}


# Computing bias ---------------------------------------------------------------
bias_LBFGSB <- do.call(rbind, Map(bias_calc, ans_LBFGSB[1:i], par[1:i], dat[1:i], 1:i))
bias_LBFGSB_priors <- do.call(rbind, Map(bias_calc, ans_LBFGSB_priors[1:i], par[1:i], dat[1:i], 1:i))
bias_ABC    <- do.call(rbind, Map(bias_calc, ans_ABC[1:i], par[1:i], dat[1:i], 1:i))
bias_ABC_priors    <- do.call(rbind, Map(bias_calc, ans_ABC_priors[1:i], par[1:i], dat[1:i], 1:i))
bias_MCMC   <- do.call(rbind, Map(bias_calc, ans_MCMC[1:i], par[1:i], dat[1:i], 1:i))

# Checking solved solutions

bias <- rbind(
  data.frame(Method = "MLE w/ L-BFGS-B", bias_LBFGSB),
  # data.frame(Method = "MLE+priors w/ L-BFGS-B", bias_LBFGSB_priors),
  data.frame(Method = "MCMC", bias_MCMC)
  # data.frame(Method = "MLE w/ ABC", bias_ABC),
  # data.frame(Method = "MLE+priors w/ ABC", bias_ABC_priors)
)

# Reshaping 
bias <- do.call(rbind, lapply(colnames(bias)[-(1:5)], function(x) {
  ans <- data.frame(parameter = x, bias[,c("Method", "index", "TreeSize", "NLeafs", "Missings",x),drop=FALSE])
  colnames(ans)[7] <- "bias"
  ans
}
))

# Dropping missings
bias <- bias[complete.cases(bias),]

# Categorial variables ---------------------------------------------------------

# Missings
bias$miss_tag <- interval_tags(bias$Missings, 1:5/ 10)

# Tree size
bias$size_tag <- interval_tags(bias$TreeSize, quantile(bias$TreeSize, na.rm = TRUE))

# NLeafs/TreeSize
bias$PropLeafs <- with(bias, NLeafs/TreeSize)
bias$PropLeafs_tag <- interval_tags(bias$PropLeafs, quantile(bias$PropLeafs, na.rm=TRUE))

#+ plotting, echo=TRUE
# Plot -------------------------------------------------------------------------
graphics.off()
sizelvls <- levels(bias$size_tag)
for (i in 1:4) {
  # Clearing plot space and creating the pdf
  
  # pdf(sprintf("playground/simulations/bias_trees_of_size_%s.pdf", sizelvls[i]))
  
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
  
  # dev.off()
  
}
