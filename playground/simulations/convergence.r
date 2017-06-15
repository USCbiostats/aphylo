#' ---
#' title: "Convergence"
#' author: "George G Vega Yon"
#' date: "`r paste('This version:', Sys.time())`"
#' ---
#' 

#+ setup,
knitr::opts_chunk$set(echo = FALSE)

#+ data-loading, cache=TRUE
load("playground/simulations/simdata3.rda")

# Tabulating error type
tabulate_error <- function(object) {
  nvalid <- length(which(sapply(object, length) > 0))
  object <- object[1:nvalid]
  
  ans_conv <- sapply(object, function(x) {
    if (inherits(x, "error")) return("Error")
    
    ans <- x[["varcovar"]]
    
    if (any(diag(ans) < 0)) "No convergence"
    else "Convergence"
  })
  
  ans_msg <- sapply(object, function(x) {
    if (inherits(x, "error")) return("Error")
    
    if (!is.na(x[["message"]])) x[["message"]]
    else "No message"
  })
  
  print(addmargins(prop.table(table(ans_msg, ans_conv)))*100, digits=2)
  
}

#+ echo=TRUE
tabulate_error(ans_MLE)
tabulate_error(ans_MAP)
tabulate_error(ans_MCMC_right_prior)
tabulate_error(ans_MCMC_wrong_prior)

gel <- lapply(ans_MCMC_right_prior, "[[", "gelman.mpsrf")


