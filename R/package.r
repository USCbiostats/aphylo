#' @importFrom Rcpp evalCpp
#' @importFrom stringr str_match_all
#' @importFrom stats var coef vcov optim optimHess
#' @importFrom utils head
#' @importFrom graphics plot contour persp legend mtext plot.new plot.window par
#'  points
#' @importFrom grDevices colors trans3d adjustcolor
#' @importFrom coda mcmc mcmc.list
#' @importFrom amcmc MCMC
#' @importFrom ape as.phylo plot.phylo
#' @importFrom ggtree ggtree geom_tiplab gheatmap
#' @importFrom ggplot2 scale_fill_manual
NULL

#' @useDynLib aphylo
NULL

# Importing from the Matrix pkg ------------------------------------------------

#' @import methods
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom Matrix Matrix rowSums colSums
#' @importMethodsFrom Matrix t
NULL

# Data documentation -----------------------------------------------------------

#' Phylogenetic Tree
#' 
#' A dataset containing the parent-offspring relations between genes.
#' 
#' @format A data frame with 326 rows and 3 variables:
#' \describe{
#'   \item{NodeId}{Integer, ID of the offspring.}
#'   \item{TypeId}{Integer, Type of node (unusued).}
#'   \item{ParentId}{Integer, ID of the parent.}
#' }
#' 
#' @source BiostatsUSC
#' @name tree
NULL

#' Experimental Data
#' 
#' A dataset containing 3 functional state of the leaf nodes. Each
#' function can have either 0 (unactive), 1 (active) or 9 (n/a).
#' 
#' @format A data frame with 365 rows and 4 variables:
#' \describe{
#'   \item{f01}{State of function 1.}
#'   \item{f02}{State of function 1.}
#'   \item{f03}{State of function 1.}
#'   \item{LeafId}{Integer, ID of the leaf.}
#' }
#' 
#' @source BiostatsUSC
#' @name experiment
NULL

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
#' @docType package
#' @name aphylo
#' 
NULL