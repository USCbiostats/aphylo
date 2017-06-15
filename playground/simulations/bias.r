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
load("playground/simulations/data_and_functions.rda")

library(ggplot2)
library(magrittr)

# Function to measure bias
bias_calci <- function(x, par0, tree) {

  if (!length(x)) return(NULL)

  # Names of the objects that will be stored
  cnames <- c("psi0", "psi1", "mu0", "mu1", "Pi")
  cnames <- c(
    "TreeSize",
    "NLeafs",
    "Missings",
    cnames,
    paste(cnames, "bias", sep="_"),
    paste(cnames, "var", sep="_")
  )
  
  # Checking if it was able to solve it or not
  if (inherits(x, "error"))
    return(structure(rep(NA, length(cnames)), names = cnames))

  # Number of offspring and internal nodes
  treesize <- length(tree$noffspring)
  nleafs   <- sum(tree$noffspring > 0)


  structure(
    c(
      treesize,
      nleafs,
      par0[6],
      x$par,
      x$par - par0[1:5],
      diag(x$varcovar)
      ),
    names = cnames
    )

}

# Function to try to load a file -ntries- times waiting -wait- seconds
# between each try.
tryLoad <- function(fn, envir = .GlobalEnv, ..., ntries = 5, wait = 120) {
  
  message("Loading file ",fn, "...", appendLF = FALSE)
  i <- 1
  while (i < ntries) {
    ans <- tryCatch(load(fn, envir = envir), error = function(e) e)
  
    if (inherits(ans, "error"))  {
      i <- i + 1
      print(ans)
      message("Error! trying again in ", wait, "seconds (",i,"/",ntries,")...")
      Sys.sleep(wait)
      next
    }
    
    break
  }
  
  message("done.")
  
}

bias_calc <- function(fn, objname) {
  env <- new.env()

  # Trying to load the file
  tryLoad(fn, envir = env)
  ids <- which(sapply(env[[objname]], length) > 0)

  message("\tComputing bias ...", appendLF = FALSE)
  ans <- Map(bias_calci, env[[objname]][ids], parameters[ids], dat[ids])
  ans <- do.call(rbind, ans)
  message("done.")

  ans <- cbind(index = ids, ans)

  ans[complete.cases(ans),,drop=FALSE]
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
bias_MLE <- bias_calc("playground/simulations/mle_estimates.rda", "ans_MLE")
bias_MAP <- bias_calc("playground/simulations/map_estimates.rda", "ans_MAP")
bias_MAP_wrong <- bias_calc("playground/simulations/map_wrong_prior_estimates.rda", "ans_MAP_wrong_prior")
bias_MCMC_right <- bias_calc("playground/simulations/mcmc_right_prior_estimates.rda", "ans_MCMC_right_prior")
bias_MCMC_wrong <- bias_calc("playground/simulations/mcmc_wrong_prior_estimates.rda", "ans_MCMC_wrong_prior")

# Checking solved solutions
common_solutions <- intersect(bias_MLE[,"index"], bias_MCMC_right[,"index"])
common_solutions <- intersect(common_solutions, bias_MAP_wrong[,"index"])
common_solutions <- intersect(common_solutions, bias_MCMC_right[,"index"])
common_solutions <- intersect(common_solutions, bias_MCMC_wrong[,"index"])

bias <- rbind(
  data.frame(Method = "MLE", bias_MLE),
  data.frame(Method = "MAP", bias_MAP),
  data.frame(Method = "MAP wrong", bias_MAP_wrong),
  data.frame(Method = "MCMC wrong", bias_MCMC_wrong),
  data.frame(Method = "MCMC", bias_MCMC_right)
)

# Categorial variables ---------------------------------------------------------

# Missings
bias$miss_tag <- interval_tags(bias$Missings, 1:5/ 10)

# Tree size
bias$size_tag <- interval_tags(bias$TreeSize, quantile(bias$TreeSize, na.rm = TRUE))

# NLeafs/TreeSize
bias$PropLeafs <- with(bias, NLeafs/TreeSize)
bias$PropLeafs_tag <- interval_tags(bias$PropLeafs, quantile(bias$PropLeafs, na.rm=TRUE))

save(bias, file = "playground/simulations/bias.rda")


