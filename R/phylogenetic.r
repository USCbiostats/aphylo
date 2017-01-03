#' @importFrom Rcpp evalCpp
#' @importFrom igraph plot.igraph layout_with_sugiyama graph_from_edgelist
#' @importFrom graphics plot
NULL

#' @useDynLib phylogenetic
NULL

# Importing from the Matrix pkg ------------------------------------------------

#' @import methods
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom Matrix Matrix rowSums colSums
#' @importMethodsFrom Matrix t
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