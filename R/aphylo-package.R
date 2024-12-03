#' @importFrom Rcpp evalCpp
#' @importFrom stats var coef vcov optim optimHess complete.cases dbeta runif
#'  as.formula update.formula predict window
#' @importFrom utils head tail getFromNamespace combn
#' @importFrom graphics plot contour persp legend mtext plot.new plot.window par
#'  points segments polygon rect yinch image abline
#' @importFrom grDevices colors trans3d adjustcolor rgb
#' @importFrom coda mcmc mcmc.list
#' @importFrom fmcmc MCMC
#' @importFrom MASS ginv
#' @importFrom xml2 as_list read_html
NULL

#' @useDynLib aphylo, .registration=TRUE
NULL

# Importing from the Matrix pkg ------------------------------------------------

#' @import methods
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom Matrix Matrix rowSums colSums
#' @importMethodsFrom Matrix t
NULL

# Data documentation -----------------------------------------------------------

#' Fake Phylogenetic Tree
#' 
#' A fake dataset containing the parent-offspring relations between genes.
#' This dataset is inteded for testing only.
#' 
#' @format A data frame with 6 rows and 2 variables:
#' \describe{
#'   \item{NodeId}{Integer, ID of the offspring.}
#'   \item{ParentId}{Integer, ID of the parent.}
#' }
#' 
#' @source BiostatsUSC
#' @name faketree
NULL

#' Fake Experimental Data
#' 
#' A fake dataset containing 2 functional state of the leaf nodes. Each
#' function can have either 0 (unactive), 1 (active) or 9 (n/a).
#' This dataset is inteded for testing only.
#' 
#' @format A data frame with 4 rows and 3 variables:
#' \describe{
#'   \item{f1}{State of function 1.}
#'   \item{f2}{State of function 1.}
#'   \item{LeafId}{Integer, ID of the leaf.}
#' }
#' 
#' @source BiostatsUSC
#' @name fakeexperiment
NULL


# Importing from R CORE packages -----------------------------------------------

release_questions <- function() {
  c(
    "Have you updated the inst/NEWS file?",
    "Have you changed the version+dates in DESCRIPTION and NEWS.md?",
    "Have you added all new files to GIT?",
    "Have you clean the vignettes file (source)?"
  )
}

#' Statistical Inference in Annotated Phylogenetic Trees
#' 
#' Statistical Inference in Annotated Phylogenetic Trees
#' 
#' @name aphylo-package
#' 
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
  options(
    aphylo_informative = FALSE,
    aphylo_reduce_pseq = TRUE
    )
}
